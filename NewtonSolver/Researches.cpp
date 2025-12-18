// Researches.cpp


#include "Researches.h"
#include <string>
#include <algorithm>

const double EPS = 1e-8;

struct ICurve2D
{
    virtual ~ICurve2D() = default;

    // значение F(x, y)
    virtual double value(double x, double y) const = 0;

    // частные производные
    virtual double d_dx(double x, double y) const = 0;
    virtual double d_dy(double x, double y) const = 0;
};
// Окружность: (x - cx)^2 + (y - cy)^2 - r^2 = 0
struct Circle : public ICurve2D
{
public: 
    double cx, cy, r;

    Circle(double cx, double cy, double r) : cx(cx), cy(cy), r(r) {}

    double value(double x, double y) const override
    {
        double dx = x - cx;
        double dy = y - cy;
        return dx * dx + dy * dy - r * r;
    }

    double d_dx(double x, double y) const override
    {
        (void)y;
        return 2.0 * (x - cx);
    }

    double d_dy(double x, double y) const override
    {
        (void)x;
        return 2.0 * (y - cy);
    }
};
// Прямая: ty = k * x + b  ->  ty - kx - b = 0
struct Line : public ICurve2D
{
public:
    double t;
    double k;
    double b;

    Line(double t, double k, double b) : t(t), k(k), b(b) {}

    double value(double x, double y) const override
    {
        return t * y - k * x - b;
    }

    double d_dx(double x, double y) const override
    {
        (void)x; (void)y;
        return -k;
    }

    double d_dy(double x, double y) const override
    {
        (void)x; (void)y;
        return t;
    }
};
// Синусоида: y = A * sin(B * x) + C  →  y - (A sin(Bx) + C) = 0
struct SineCurve : public ICurve2D
{
public:
    double A, B, C;

    SineCurve(double A, double B, double C) : A(A), B(B), C(C) {}

    double value(double x, double y) const override
    {
        return y - (A * sin(B * x) + C);
    }

    double d_dx(double x, double y) const override
    {
        (void)y;
        // d/dx [y - (A sin(Bx) + C)] = -A B cos(Bx)
        return -A * B * cos(B * x);
    }

    double d_dy(double x, double y) const override
    {
        (void)x; (void)y;
        // d/dy [y - ...] = 1
        return 1.0;
    }
};
class System2D
{
private:
    vector<shared_ptr<ICurve2D>> equations_;
    vector<double> weights_; 
    vector<NewtonSolver::Equation> newtonEquations_;
public:
    // Добавить произвольную кривую
    System2D& addCurve(shared_ptr<ICurve2D> curve)
    {
        equations_.push_back(move(curve));
        weights_.push_back(1.0);
        return *this;
    }

    System2D& addCircle(double cx, double cy, double r)
    {
        return addCurve(make_shared<Circle>(cx, cy, r));
    }

    System2D& addLine(double t, double k, double b)
    {
        return addCurve(make_shared<Line>(t, k, b));
    }

    System2D& addLine(double k, double b)
    {
        return addLine(1.0, k, b);
    }
    System2D& addPerpendicularToOx(double x0)
    {
        return addLine(0.0, 1.0, -x0);
    }
    System2D& addSine(double A, double B, double C)
    {
        return addCurve(std::make_shared<SineCurve>(A, B, C));
    }
    int equationsCount() const { return static_cast<int>(equations_.size()); }

    void setWeight(int index, double w)
    {
        if (index >= 0 && index < (int)weights_.size())
            weights_[index] = w;
    }

    void setWeights(const double* w, int count)
    {
        for (int i = 0; i < count && i < (int)weights_.size(); ++i)
            weights_[i] = w[i];
    }

    NewtonSolver::Equation* equationsFor(NewtonSolver& solver)
    {
        (void)solver;

        newtonEquations_.resize(equations_.size());

        for (size_t i = 0; i < equations_.size(); ++i)
        {
            auto curve = equations_[i];
            double w = weights_.empty() ? 1.0 : weights_[i];

            // F(x) = w * F_curve(x)
            newtonEquations_[i].F =
                [curve, w](const double* x, int n) -> double
                {
                    (void)n;
                    return w * curve->value(x[0], x[1]);
                };

            // dF/dx_j = w * dF_curve/dx_j
            newtonEquations_[i].dF =
                [curve, w](const double* x, int n, int j) -> double
                {
                    (void)n;
                    if (j == 0) return w * curve->d_dx(x[0], x[1]);
                    if (j == 1) return w * curve->d_dy(x[0], x[1]);
                    return 0.0;
                };
        }

        return newtonEquations_.data();
    }
};
void PrintHeader(const char* title)
{
    cout << string(50, '=') << endl;
    cout << title << endl;
    cout << string(50, '=') << endl;
}

void RunForStarts(NewtonSolver& solver,
    const double starts[][2],
    int count,
    bool useAnalyticJacobian,
    bool useVariant1,
    const char* title)
{
    int s;

    PrintHeader(title);

    string cleanTitle = title;
    auto is_bad_char = [](char c) {
        return c == ' ' || c == ':' || c == ',' || c == '(' || c == ')' || c == '/';
        };
    replace_if(cleanTitle.begin(), cleanTitle.end(), is_bad_char, '_');

    string::iterator new_end = unique(cleanTitle.begin(), cleanTitle.end(),
        [](char a, char b) { return a == '_' && b == '_'; });
    cleanTitle.erase(new_end, cleanTitle.end());

    for (s = 0; s < count; ++s)
    {
        // задаём начальное приближение
        solver.SetX(0, starts[s][0]);
        solver.SetX(1, starts[s][1]);
        string filename = "result_" + cleanTitle + "_start" + to_string(s) + ".csv";
        cout << "Старт #" << s << " -> " << filename << endl;
        // решаем
        solver.NewtonSolve(useAnalyticJacobian, useVariant1, filename);

        cout << "Начальное приближение: ("
            << starts[s][0] << ", " << starts[s][1] << ")"
            << " - точка пересечения: (" << scientific << setprecision(15)
            << solver.GetX(0) << ", " << solver.GetX(1) << ")"
            << endl << endl << fixed << setprecision(1);
    }

    cout << std::string(50, '=') << std::endl << std::endl;
}

void CirclePlusLineResearch()
{
    double x1 = 0.0, y1 = 0.0, r1 = 1.0;
    double x2 = 3, y2 = 0.0, r2 = 1.0;
    double k = 2, b = -0.5;

    NewtonSolver solver;
    solver.AllocateMemory(2, 3);
    solver.SetEps1(1e-15);
    solver.SetEps2(1e-15);
    solver.SetMaxIter(50);

    System2D system;
    system.addCircle(x1, y1, r1)
          .addCircle(x2, y2, r2)
          //.addPerpendicularToOx(0.5);
        ;
    solver.SetSystem(system.equationsCount(),
        system.equationsFor(solver));

    double starts[4][2] = {
        {0.5, 0.5},
        {1.5, 0.0},
        {1.5, 1.0},
        {0.0, 0.0}
    };
    
    double weightSets[][3] = {
        {1.0, 1.0, 1.0},   // без взвешивания
        //{1.0, 1.0, 1.0},  // первая окружность важнее
        //{1.0, 1.0, 1.0},  // вторая окружность важнее
        //{1.0, 1.0, 1.0}   // прямая важнее
    };
    const int numWeightSets = sizeof(weightSets) / sizeof(weightSets[0]);

    for (int w = 0; w < numWeightSets; ++w)
    {
        // задаём веса для трёх уравнений
        system.setWeights(weightSets[w], 3);

        cout << "Веса уравнений: ["
            << weightSets[w][0] << ", "
            << weightSets[w][1] << ", "
            << weightSets[w][2] << "]" << endl;

        // пересобираем систему с учётом весов
        solver.SetSystem(system.equationsCount(),
            system.equationsFor(solver));

        RunForStarts(solver, starts, 4,
            true, true,
            "Вариант 1: аналитический Якобиан");

        RunForStarts(solver, starts, 4,
            true, false,
            "Вариант 2: аналитический Якобиан");

        RunForStarts(solver, starts, 4,
            false, false,
            "Вариант 6: численный Якобиан");
    }
}
void ThreeLineWithWeight()
{
    // Три попарно пересекающиеся прямые:
     // L1: y = -2x + 4
     // L2: y = x + 1
     // L3: y = 1        (горизонтальная)
    double k[3] = { -2.0,  1.0, 0.0 };
    double b[3] = { 4.0,  1.0, 1.0 };

    double starts[4][2] = {
        {0.0, 0.0},
        {1.0, 1.0},
        {2.0, 1.5},
        {3.0, 2.0}
    };
    // Наборы весов для пары прямых: [w1, w2]
    double weightSets[][2] = {
        {1.0, 1.0},       // без взвешивания
        {10000.0, 1.0},   // первая прямая важнее
        {1.0, 10000.0}    // вторая прямая важнее
    };
    const int numWeightSets = sizeof(weightSets) / sizeof(weightSets[0]);

    NewtonSolver solver;
    solver.AllocateMemory(2, 2);   // 2 переменные (x, y), 2 уравнения (2 прямые)
    solver.SetEps1(1e-15);
    solver.SetEps2(1e-15);
    solver.SetMaxIter(50);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = i + 1; j < 3; ++j)
        {
            cout << std::string(80, '*') << endl;
            cout << "Пара прямых L" << (i + 1) << " и L" << (j + 1) << endl;
            cout << "L" << (i + 1) << ": y = " << k[i] << " * x + " << b[i] << endl;
            cout << "L" << (j + 1) << ": y = " << k[j] << " * x + " << b[j] << endl;
            cout << std::string(80, '*') << endl;

            for (int w = 0; w < numWeightSets; ++w)
            {
                System2D system;
                system.addLine(k[i], b[i])   // уравнение 0
                    .addLine(k[j], b[j]);  // уравнение 1

                // задаём веса для двух уравнений
                system.setWeights(weightSets[w], 2);

                cout << "Веса уравнений: ["
                    << weightSets[w][0] << ", "
                    << weightSets[w][1] << "]" << endl;

                // пересобираем систему с учётом весов
                solver.SetSystem(system.equationsCount(),
                    system.equationsFor(solver));
                
                // Вариант 1: аналитический Якобиан
                RunForStarts(solver, starts, 4,
                    true, true,
                    "Вариант 1: аналитический Якобиан");
                
                // Вариант 2: аналитический Якобиан
                RunForStarts(solver, starts, 4,
                    true, false,
                    "Вариант 2: аналитический Якобиан");

                // Вариант 6: численный Якобиан
                RunForStarts(solver, starts, 4,
                    false, false,
                    "Вариант 6: численный Якобиан");

                cout << std::string(80, '#') << endl << endl;
            }
        }
    }
}
void SinusoidLineResearch()
{
    // Параметры синусоиды и прямой
    // y = sin(x)
    double A = 1.0, B = 1.0, C = 0.0;

    // y = kx + b (почти горизонтальная, чтобы было несколько пересечений)
    double k = 0.02;
    double b = 0.0;

    NewtonSolver solver;
    solver.AllocateMemory(2, 2);   // 2 переменные (x, y), 2 уравнения
    solver.SetEps1(1e-15);
    solver.SetEps2(1e-15);
    solver.SetMaxIter(50);

    System2D system;
    system.addSine(A, B, C)   // F1: y - sin(x) = 0
        .addLine(k, b);     // F2: y - (k x + b) = 0

    solver.SetSystem(system.equationsCount(),
        system.equationsFor(solver));

    // Набор начальных приближений:
    // берем точки слева, около разных пересечений и максимумов/минимумов синуса
    double starts[8][2] = {
        { 0.01, 0.0 },
        { 2.0, 0.02 },
        { 0.5,  0.5},
        { 4.0,  0.0},
        { 6.0,  0.0},
        { 8.0,  0.0}
    };

    // Можно исследовать, например, вариант 2 и 6,
    // как в других исследованиях
    RunForStarts(solver, starts, 6,
        true, true,
        "Синусоида + прямая, вариант 1 (аналитический Якобиан)");

    RunForStarts(solver, starts, 6,
        true, false,
        "Синусоида + прямая, вариант 2 (аналитический Якобиан)");

    RunForStarts(solver, starts, 6,
        false, false,
        "Синусоида + прямая, вариант 6 (численный Якобиан)");
}