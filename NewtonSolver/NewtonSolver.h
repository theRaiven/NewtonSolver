#pragma once
#include <iostream>
#include <functional>
#include <fstream>
using namespace std;

using real = double;
using realS = double;
const int accuracy = 15;

template<typename T>
constexpr T relativeEPS();

template<>
constexpr float relativeEPS<float>() { return 1e-34f; }

template<>
constexpr double relativeEPS<double>() { return 1e-300; }

void SolveGauss(real** A, real* b, real*& x, int& n);
using Func = function<double(const double*, int)>;
using DerivativeFunc = double(*)(const double* x, int n, int j);


class NewtonSolver
{
private:
    int n;          // число неизвестных
    int m;          // число уравнений

    real eps1;      
    real eps2;      
    int maxIter;    

    real* x;        // текущее приближение x^k
    real* Fx;       // значение F(x^k)
    real* deltaX;   // Δx^k

    real** J;       // матрица Якоби (m × n)
    
    // служебные буферы для квадратной СЛАУ
    real** Aquad; // квадратная матрица
    real* bquad;  // квадратный вектор RHS
    real* xquad;  // решение квадратной задачи
    
    int* selectedVars;       // размер n
    int* selectedEquations;  // размер m
public:
    struct Equation 
    {
        Func F;  // F_i(x)
        function<double(const double*, int, int)> dF; // ∂F_i/∂x_j
    };
private:
    Equation* eq;

public:
    NewtonSolver();
    ~NewtonSolver();

	// методы доступа к параметрам 
    void SetX(int i, double val) { x[i] = val; }
    double GetX(int i) const { return x[i]; }
    void SetEps1(double val) { eps1 = val; }
    void SetEps2(double val) { eps2 = val; }
    void SetMaxIter(int val) { maxIter = val; }
    int GetN() const { return n; }
    void AllocateMemory(int numVars, int numEqs)
    {
        n = numVars;
        m = numEqs;

        // Векторы
        x = new real[n]();         // текущее приближение
        deltaX = new real[n]();    // Δx^k
        Fx = new real[m]();        // F(x^k)

        // Якобиан m x n
        J = new real * [m];
        for (int i = 0; i < m; i++)
            J[i] = new real[n]();

        // Квадратная СЛАУ n x n
        Aquad = new real * [n];
        for (int i = 0; i < n; i++)
            Aquad[i] = new real[n]();

        bquad = new real[n]();     // RHS
        xquad = new real[n]();     // решение квадратной задачи

        // Выбор переменных и уравнений
        selectedVars = new int[n];
        selectedEquations = new int[m];

        // Уравнения еще не зарегистрированы
        eq = nullptr;
    }

    // служеные методы 
    bool LoadConfig(const string configFile, const string startFile);
    void SaveResult(const string fileName);
    void PrintState();
    double Norm(real* a, int n);
    
    // регистрация системы уравнений
    void SetSystem(int mEquations, Equation* equations);
    void ComputeF();      
    bool IsJacobianSingular(double tol = 1e-15);
    void ComputeJacobian();
    void ComputeJacobianNumeric(real h = 1e-8);

    // решатели 
    void SelectVariablesVariant1(int* selectedVars);
    void SelectEquationsVariant2(int* selectedEquations);

    void FormSquareMatrixVariant1(const int* selectedVars);
    void FormSquareMatrixVariant2(const int* selectedEquations);
    
    void SolveDeltaX(const int* selected, bool isVariant1);
    double FindBeta();
    void NewtonSolve(bool useAnalyticJacobian = true, bool selectVars = true); // и собственно, тот, из-за кого мы все здесь собрались 

};