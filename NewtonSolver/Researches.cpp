#include "Researches.h"

const double EPS = 1e-8;

// Евклидово расстояние между точками
double Dist(const double* a, const double* b, int n = 2)
{
    double sum = 0;
    for (int i = 0; i < n; i++) sum += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(sum);
}

// Проверка, есть ли точка в списке (с учётом EPS)
bool IsNewSolution(const vector<vector<double>>& sols, const double* x)
{
    for (const auto& s : sols)
        if (Dist(s.data(), x) < EPS) return false;
    return true;
}

// Создание системы для двух окружностей
NewtonSolver::Equation* CreateCircleSystem(double x1, double y1, double r1,
    double x2, double y2, double r2)
{
    auto eq = new NewtonSolver::Equation[2];

    eq[0].F = [=](const double* x, int) -> double {
        return (x[0] - x1) * (x[0] - x1) + (x[1] - y1) * (x[1] - y1) - r1 * r1;
        };
    eq[1].F = [=](const double* x, int) -> double {
        return (x[0] - x2) * (x[0] - x2) + (x[1] - y2) * (x[1] - y2) - r2 * r2;
        };

    eq[0].dF = [=](const double* x, int, int j) -> double {
        return j == 0 ? 2 * (x[0] - x1) : 2 * (x[1] - y1);
        };
    eq[1].dF = [=](const double* x, int, int j) -> double {
        return j == 0 ? 2 * (x[0] - x2) : 2 * (x[1] - y2);
        };

    return eq;
}
NewtonSolver::Equation* CreateCirclePlusLineSystem(double x1, double y1, double r1,
    double x2, double y2, double r2,
    double k, double b) // параметры прямой y = k*x + b
{
    auto eq = new NewtonSolver::Equation[3];

    // Первая окружность
    eq[0].F = [=](const double* x, int) -> double {
        return (x[0] - x1) * (x[0] - x1) + (x[1] - y1) * (x[1] - y1) - r1 * r1;
        };
    eq[0].dF = [=](const double* x, int, int j) -> double {
        return j == 0 ? 2 * (x[0] - x1) : 2 * (x[1] - y1);
        };

    // Вторая окружность
    eq[1].F = [=](const double* x, int) -> double {
        return (x[0] - x2) * (x[0] - x2) + (x[1] - y2) * (x[1] - y2) - r2 * r2;
        };
    eq[1].dF = [=](const double* x, int, int j) -> double {
        return j == 0 ? 2 * (x[0] - x2) : 2 * (x[1] - y2);
        };

    // Прямая y = k*x + b
    eq[2].F = [=](const double* x, int) -> double {
        return x[1] - k * x[0] - b;
        };
    eq[2].dF = [=](const double* x, int, int j) -> double {
        return j == 0 ? -k : 1;
        };

    return eq;
}

// Основной тест
void CircleIntersectionResearch()
{
    // Пример: две окружности
    double x1 = 0, y1 = 0, r1 = 1;
    double x2 = 2.0, y2 = 0, r2 = 1; // не пересекаются
    double k = 0.5, b = -0.5;         // прямая y = k*x + b

    NewtonSolver solver;
    auto eq = CreateCirclePlusLineSystem(x1, y1, r1, x2, y2, r2, k, b);
    solver.AllocateMemory(2, 3);   // 2 переменные, 3 уравнения
    solver.SetSystem(3, eq);

    vector<vector<double>> initialGuesses1 = {
         {0.5, 0.5, 0.0},   // не на осях симметрии
    {1.0, 0.0, 0.0},   // на линии центров
    {1.0, 0.5, 0.0},   // на перпендикулярной оси на равном расстоянии
    {2.0, 0.0, 0.0}    // внутри первой окружности
    };
    vector<vector<double>> initialGuesses = {
       {0.5, 0.5},
       {1.0, 0.0},
       {1.0, 0.5},
       {2.0, 0.0}
    };
    vector<vector<double>> solutions;

    solver.SetEps1(1e-15);
    solver.SetEps2(1e-15);
    solver.SetMaxIter(50);

    int n = solver.GetN();

    for (auto& start : initialGuesses)
    {
        for (int i = 0; i < n; i++) solver.SetX(i, start[i]);

        cout << "Начальное приближение: ("
            << start[0] << ", " << start[1] << ")\n";

        solver.NewtonSolve(false, true); // анал, вар2

        // формируем массив x для проверки
        double xArr[2] = { solver.GetX(0), solver.GetX(1) };

        IsNewSolution(solutions, xArr);
        solutions.push_back({ xArr[0], xArr[1] });
        cout << "Найдена новая точка пересечения: " << fixed << setprecision(8) << xArr[0] << ", " << xArr[1] << endl;
        cout << string(50, '-') << endl << endl;
    }

    if (solutions.empty())
        cout << "Решений нет.\n";
    else
        cout << "Всего найдено точек пересечения: " << solutions.size() << endl;

    delete[] eq;
}