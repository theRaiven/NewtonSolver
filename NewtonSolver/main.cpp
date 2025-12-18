#include "NewtonSolver.h"
#include "Researches.h"



int main()
{
	setlocale(LC_ALL, "Russian");

    cout << "=== ЗАПУСК ИССЛЕДОВАНИЙ ===" << endl;

    // ---------------------------------------------------------
    // 1. ИССЛЕДОВАНИЕ ОКРУЖНОСТЕЙ (Часть А, Б, В)
    // ---------------------------------------------------------

    // а) Не пересекаются (Центры разнесены на 3, радиусы по 1)
    cout << "\n>> Окружности: Не пересекаются..." << endl;
    CirclePlusLineResearch(0.0, 0.0, 1.0, 3.0, 0.0, 1.0);

    // б) Касаются (одна точка) (Центры разнесены на 2)
    cout << "\n>> Окружности: Касание..." << endl;
    CirclePlusLineResearch(0.0, 0.0, 1.0, 2.0, 0.0, 1.0);

    // в) Пересекаются (две точки) (Центры разнесены на 1)
    cout << "\n>> Окружности: Пересечение (2 точки)..." << endl;
    CirclePlusLineResearch(0.0, 0.0, 1.0, 1.0, 0.0, 1.0);

    cout << "\n>> Окружности + Прямая..." << endl;
    CirclePlusLineResearch(0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.5);
    CirclePlusLineResearch(0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.5, 0.5);
    ResearchVerticalLine();
    cout << "\n>> Три прямые..." << endl;

    double k_dummy[3] = { -2.0, 1.0, 0.0 };
    double b_dummy[3] = { 4.0, 1.0, 1.0 };

    double weights[][2] = {
        {1.0, 1.0},      // Равные
        {10000.0, 1.0},  // Первая важнее
        {1.0, 10000.0}   // Вторая важнее
    };

    double* weights_ptrs[3] = { weights[0], weights[1], weights[2] };
    ThreeLineWithWeight(k_dummy, b_dummy, weights_ptrs);

    cout << "\n>> Синусоида..." << endl;
    SinusoidLineResearch();

    cout << "\nВсе готово." << endl;
    return 0;
}

