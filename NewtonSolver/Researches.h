// Researches.h


#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include "NewtonSolver.h"


void CirclePlusLineResearch(double x1, double y1, double r1, double x2, double y2, double r2, double k, double b);
void CirclePlusLineResearch(double x1, double y1, double r1, double x2, double y2, double r2);
void ThreeLineWithWeight(double* k, double* b, double** weightSets);
void SinusoidLineResearch();
void ResearchVerticalLine();