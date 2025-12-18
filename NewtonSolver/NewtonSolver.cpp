// NewtonSolver.cpp 
//

#include "NewtonSolver.h"
#include <iomanip>

void SolveGauss(real** A, real* b, real*& x, int& n)
{
    try
    {
        const real EPS = relativeEPS<real>();

        for (int i = 0; i < n; i++)
        {
            int i0 = i;
            real maxElem = fabs(A[i][i]);

            for (int j = i + 1; j < n; j++)
            {
                if (fabs(A[j][i]) > maxElem)
                {
                    maxElem = fabs(A[j][i]); 
                    i0 = j;
                }
            }

            swap(A[i], A[i0]); 
            swap(b[i], b[i0]);

            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("Система вырождена");
            }


            for (int j = i + 1; j < n; j++)
            {
                real m{ A[j][i] / A[i][i] }; 
                for (int k = i; k < n; k++)
                {
                    A[j][k] -= m * A[i][k]; 
                } 
                b[j] -= m * b[i]; 
            }
        }


        for (int i = n - 1; i >= 0; i--)
        {
            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("Система вырождена");
            }
            realS sum{ 0.0 };
            for (int j = i + 1; j < n; j++)
            {
                sum += A[i][j] * x[j]; 
            }
            x[i] = (b[i] - sum) / A[i][i]; 
        } 
    }
    catch (std::exception& e)
    {
        std::cout << "Ошибка! " << e.what() << std::endl;
    }
}

NewtonSolver::NewtonSolver() 
{
    n = m = 0;
    eps1 = eps2 = 0.0;
    maxIter = 0;

    x = nullptr;
    Fx = nullptr;
    deltaX = nullptr;

    eq = nullptr;

    J = nullptr;
    Aquad = nullptr;
    bquad = nullptr;
    xquad = nullptr;
}
NewtonSolver::~NewtonSolver()
{
    if (x) delete[] x;
    if (Fx) delete[] Fx;
    if (deltaX) delete[] deltaX;

    if (eq) delete[] eq;

    if (J)
    {
        for (int i = 0; i < m; i++)
            delete[] J[i];
        delete[] J;
    }

    if (Aquad)
    {
        for (int i = 0; i < n; i++)
            delete[] Aquad[i];
        delete[] Aquad;
    }

    if (bquad) delete[] bquad;
    if (xquad) delete[] xquad;
    if (selectedVars) delete[] selectedVars;
    if (selectedEquations) delete[] selectedEquations;
}
bool NewtonSolver::LoadConfig(const string configFile, const string startFile)
{
    ifstream conf(configFile);
    if (!conf.is_open())
    {
        cout << "Ошибка: не удалось открыть файл конфигурации\n";
        return false;
    }

    conf >> eps1 >> eps2 >> maxIter;
    conf >> m >> n;
    conf.close();

    x = new real[n];
    deltaX = new real[n];

    if (Fx) delete[] Fx;
    Fx = new real[m];

    // матрица Якоби m×n
    J = new real * [m];
    for (int i = 0; i < m; i++)
        J[i] = new real[n];

    // квадратная СЛАУ всегда размером n×n
    Aquad = new real * [n];
    for (int i = 0; i < n; i++)
        Aquad[i] = new real[n];

    bquad = new real[n];
    xquad = new real[n];

    // загружаем начальное приближение
    ifstream start(startFile);
    if (!start.is_open())
    {
        cout << "Ошибка: не удалось открыть файл начального приближения\n";
        return false;
    }

    for (int i = 0; i < n; i++)
        start >> x[i];

    selectedVars = new int[n];
    selectedEquations = new int[m];

    start.close();
    return true;
}
void NewtonSolver::SaveResult(const string fileName)
{
    ofstream fout(fileName);
    for (int i = 0; i < n; i++)
        fout << x[i] << endl;
    fout.close();
}
void NewtonSolver::PrintState()
{
    cout << "eps1 = " << eps1 << endl;
    cout << "eps2 = " << eps2 << endl;
    cout << "maxIter = " << maxIter << endl;
    cout << "m = " << m << "  n = " << n << endl;

    cout << "x0 = ";
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << endl;
}
double NewtonSolver::Norm(real* a, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
        norm += a[i] * a[i];
    return sqrt(norm);
}

void NewtonSolver::ComputeF()
{
    for (int i = 0; i < m; i++)
        Fx[i] = eq[i].w * eq[i].F(x, n);
}
void NewtonSolver::SetSystem(int mEquations, Equation* equations)
{
    m = mEquations;

    if (eq) delete[] eq;
    eq = new Equation[m];
    for (int i = 0; i < m; i++) eq[i] = equations[i];

    if (Fx) delete[] Fx;
    Fx = new real[m];

    if (J) {
        for (int i = 0; i < m; i++) delete[] J[i];
        delete[] J;
    }

    J = new real * [m];
    for (int i = 0; i < m; i++)
        J[i] = new real[n];
}

void NewtonSolver::ComputeJacobian()
{
   // cout << endl;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            J[i][j] = eq[i].w * eq[i].dF(x, n, j);
            //cout << J[i][j] << ' ';
        }
        //cout << endl;
    }
}
void NewtonSolver::ComputeJacobianNumeric(double hStep)
{
    if (!J || !x) return;
    const double eps = std::numeric_limits<double>::epsilon();

    vector<double> Fx_old(m);
    for (int i = 0; i < m; ++i) Fx_old[i] = Fx[i];

    for (int j = 0; j < n; ++j)
    {
        double hLocal = sqrt(eps) * (1.0 + fabs(x[j]));
        if (hStep > hLocal) hLocal = hStep;

        double old = x[j];
        x[j] += hLocal;
        ComputeF(); // Fx обновляется
        for (int i = 0; i < m; ++i)
            J[i][j] = (Fx[i] - Fx_old[i]) / hLocal;
        x[j] = old;
    }
}

void NewtonSolver::SelectVariablesVariant1(int* selectedVars)
{
    if (m > n) throw runtime_error("m > n. Не делай так");

    for (int k = 0; k < m; k++)
    {
        double maxVal = -1.0;
        int maxIndex = -1;

        for (int j = 0; j < n; j++)
        {
            // проверяем, не выбрана ли j раньше
            bool alreadySelected = false;
            for (int t = 0; t < k; t++)
                if (selectedVars[t] == j) { alreadySelected = true; break; }
            if (alreadySelected) continue;

            // максимальное по строкам J[:, j]
            double colMax = 0.0;
            for (int i = 0; i < m; i++)
            {
                double val = fabs(J[i][j]);
                if (val > colMax) colMax = val;
            }

            if (colMax > maxVal)
            {
                maxVal = colMax;
                maxIndex = j;
            }
        }

        selectedVars[k] = maxIndex; // выбираем переменную
    }
}
void NewtonSolver::SelectEquationsVariant2(int* selectedEquations)
{
    if (m < n) throw runtime_error("m < n. Не делай так");
    for (int k = 0; k < n; k++)
    {
        double maxVal = -1.0;
        int maxIndex = -1;

        for (int i = 0; i < m; i++)
        {
            // проверяем, не выбрано ли i раньше
            bool alreadySelected = false;
            for (int t = 0; t < k; t++)
                if (selectedEquations[t] == i) { alreadySelected = true; break; }
            if (alreadySelected) continue;

            double val = fabs(Fx[i]); // текущая величина F_i(x)
            if (val > maxVal)
            {
                maxVal = val;
                maxIndex = i;
            }
        }

        selectedEquations[k] = maxIndex; // выбираем уравнение
    }
}

void NewtonSolver::FormSquareMatrixVariant1(const int* selectedVars)
{
    for (int i = 0; i < n; i++)
    {
        bquad[i] = -Fx[i]; // RHS = -F_i
        for (int j = 0; j < n; j++)
        {
            // если переменная j выбрана, копируем соответствующий элемент из J
            bool isSelected = false;
            for (int t = 0; t < n; t++)
            {
                if (selectedVars[t] == j) { isSelected = true; break; }
            }
            Aquad[i][j] = isSelected ? J[i][j] : 0.0;
        }
    }
}
void NewtonSolver::FormSquareMatrixVariant2(const int* selectedEquations)
{
    for (int i = 0; i < n; i++)
    {
        int eqIdx = selectedEquations[i];   // выбранное уравнение
        bquad[i] = -Fx[eqIdx];             // RHS = -F[eqIdx]
        for (int j = 0; j < n; j++)
        {
            Aquad[i][j] = J[eqIdx][j];     // ну и копируем все переменные
        }
    }
}

void NewtonSolver::SolveDeltaX(const int* selected, bool isVariant1)
{

    SolveGauss(Aquad, bquad, xquad, n);
    for (int j = 0; j < n; j++) deltaX[j] = 0.0;

    // копируем решение в deltaX
    if (isVariant1)
    {
        // вариант 1: выбранные переменные
        for (int i = 0; i < n; i++) deltaX[selected[i]] = xquad[i];
    }
    else
    {
        // вариант 2: выбранные уравнения
        for (int i = 0; i < n; i++) deltaX[i] = xquad[i]; // копируем решение
	}
}

double NewtonSolver::FindBeta()
{
    double beta = 1.0;
    double normFx = Norm(Fx, m);  // ||F^k||
    double normFv = 0.0;

    real* xTemp = new real[n];
    for (int i = 0; i < n; i++) xTemp[i] = x[i]; // сохраняем x^k

    while (true)
    {
        // пробуем шаг
        for (int i = 0; i < n; i++)
            x[i] = xTemp[i] + beta * deltaX[i];

        ComputeF();
        normFv = Norm(Fx, m);

        if (normFv <= normFx || beta < eps1) break;

        beta /= 2.0;
    }

    // если beta < eps1, предупреждение
    if (beta < eps1) {
        cout << "Предупреждение: beta < eps1, выход из поиска" << endl;
        beta = 0.0;
    }

    // восстанавливаем x^k, но шаг Δx ещё не применён
    for (int i = 0; i < n; i++)
        x[i] = xTemp[i];

    delete[] xTemp;
    return beta;
}
void NewtonSolver::NewtonSolve(bool useAnalyticJacobian, bool selectVars, string logFileName)
{
    ofstream logFile;

    if (!logFileName.empty())
    {
        logFile.open(logFileName, std::ios::app);
        logFile.imbue(std::locale("C"));
        logFile << scientific << setprecision(16);
        logFile << "Iter;Beta;NormF";
        for (int i = 0; i < n; i++)
        {
            logFile << ";x_" << i;
        }
        logFile << endl;
    }
    // начальное вычисление F(x^0)
    //cout << string(50, '=') << endl;
    try
    {
        ComputeF();
        double normF0 = Norm(Fx, m);  // норма начального F
        double normFk = normF0;

        if (logFile.is_open())
        {
            logFile << 0 << ";" << 0.0 << ";" << normF0;
            for (int i = 0; i < n; i++) logFile << ";" << x[i];
            logFile << endl;
        }

        int k;
        for (k = 0; k < maxIter; k++)
        {
            //cout << "Итерация " << k << endl;
            // 1. Вычисляем Якобиан
            if (useAnalyticJacobian)
                ComputeJacobian();       // аналитический
            else
                ComputeJacobianNumeric();

            // 2. Выбираем переменные/уравнения
            if (selectVars) SelectVariablesVariant1(selectedVars);
            else SelectEquationsVariant2(selectedEquations);
            for (int i = 0; i < n; i++) selectedVars[i] = i;
            // 3. Формируем квадратную систему
            if (selectVars) FormSquareMatrixVariant1(selectedVars);
            else FormSquareMatrixVariant2(selectedEquations);

            if (IsJacobianSingular())
            {
                for (int i = 0; i < n; i++)
                    x[i] += 1e-6;  // небольшой сдвиг
                ComputeF(); // пересчёт F после смещения
            }

            // 4. Решаем Δx^k
            SolveDeltaX(m <= n ? selectedVars : selectedEquations, m <= n);

            // 5. Поиск β
            double beta = FindBeta();
            // cout << "beta = " << beta << endl;

             // 6. Обновляем x^k+1
            for (int i = 0; i < n; i++)
                x[i] += beta * deltaX[i];

            // 7. Вычисляем F(x^k+1) и норму
            ComputeF();
            normFk = Norm(Fx, m);

            if (logFile.is_open())
            {
                logFile << (k + 1) << ";" << beta << ";" << normFk;
                for (int i = 0; i < n; i++)
                {
                    logFile << ";" << x[i];
                }
                logFile << endl;
            }

            // 8. Критерии выхода
            if (beta < eps1)
            {
                cout << "Выход: beta < eps1" << endl;
                break;
            }
            if (normFk / normF0 < eps2)
            {
                cout << "Выход: ||F^k|| / ||F^0|| < eps2" << endl;
                break;
            }
        }
        if (logFile.is_open()) logFile.close();
        if (normFk > 1e-1) // можно настраивать порог
        {
            cout << "Решений нет, ||F(x)|| = " << normFk << endl;
            for (int i = 0; i < n; i++)
            {
                //x[i] = std::numeric_limits<double>::quiet_NaN();
                // cout << "x[" << i << "] = NaN" << endl;
            }
        }
        cout << "Цикл Ньютона завершен на " << k << "-й итерации" << endl;/*
        cout << "\tAquad:";
        for (int i = 0; i < n; i++) { cout << endl << "\t\t"; for (int j = 0; j < n; j++) cout << Aquad[i][j] << " "; }
        cout << "\n\tbquad: ";
        for (int i = 0; i < n; i++) cout << bquad[i] << " "; cout << endl;

        cout << "\tНорма F = " << normFk << endl;
        cout << "\tdeltaX = ";
        for (int i = 0; i < n; i++) cout << deltaX[i] << " ";
        cout << "\n\n\tРешение x: ";
        for (int i = 0; i < n; i++)
            cout << x[i] << " ";
        cout << endl << endl;*/
    }
    catch (std::exception& e)
    {
        if (logFile.is_open()) logFile.close();
        cout << "Ошибка в методе Ньютона: " << e.what() << endl;
	}
}
bool NewtonSolver::IsJacobianSingular(double tol)
{
    for (int i = 0; i < n; i++)
        if (fabs(Aquad[i][i]) < tol) return true;
    return false;
}