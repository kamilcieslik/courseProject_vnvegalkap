//
// Created by mrfarinq on 16.06.17.
//

#include <fstream>
#include <random>
#include <climits>
#include <algorithm>
#include "TravellingSalesmanProblem.h"

TravellingSalesmanProblem::TravellingSalesmanProblem() : amountOfCities(0), arrayOfMatrixOfCities(nullptr),
                                                         optimalWay_Solution(nullptr) {
}

TravellingSalesmanProblem::~TravellingSalesmanProblem() {
    DeleteTravellingSalesman();
}

void TravellingSalesmanProblem::DeleteTravellingSalesman() {
    for (auto i = 0; i < amountOfCities; i++) {
        delete[] arrayOfMatrixOfCities[i];
    }
    delete[] arrayOfMatrixOfCities;
    arrayOfMatrixOfCities = nullptr;

    if (optimalWay_Solution != nullptr) {
        delete[] optimalWay_Solution;
        optimalWay_Solution = nullptr;
    }
}

void TravellingSalesmanProblem::LoadArrayOfMatrixOfCities(long long int **_cities, int _amountOfCities,
                                                          std::string _fileName, std::string _graphType) {
    if(arrayOfMatrixOfCities!=nullptr){
        for (int i = 0; i < amountOfCities; i++)
            delete[] arrayOfMatrixOfCities[i];
        delete[] arrayOfMatrixOfCities;
    }

    amountOfCities=_amountOfCities;

    arrayOfMatrixOfCities = new long long int *[amountOfCities];
    for (int i = 0; i < amountOfCities; i++)
        arrayOfMatrixOfCities[i] = new long long int[amountOfCities];

    for ( int i = 0; i < amountOfCities; i ++ )
    {
        for ( int j = 0; j < amountOfCities; j ++ )
        {
            arrayOfMatrixOfCities[i][j] = _cities[i][j];
        }
    }

    fileName=_fileName;
    graphType=_graphType;
}

void TravellingSalesmanProblem::GenerateRandomCities(int amountOfCities, int maxDistanceBetweenCity) {
    if (arrayOfMatrixOfCities != nullptr)
        DeleteTravellingSalesman();

    if (amountOfCities == 0) {
        std::cout << "Podaj ilość miast: ";
        std::cin >> this->amountOfCities;
        if (this->amountOfCities < 1) {
            throw std::invalid_argument("Liczba miast nie może być mniejsza od 1.");
        }

        arrayOfMatrixOfCities = new long long int *[this->amountOfCities];
        for (auto i = 0; i < this->amountOfCities; i++) {
            arrayOfMatrixOfCities[i] = new long long int[this->amountOfCities];
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist_distancesBetweenCities(1, maxDistanceBetweenCity);

        for (auto i = 0; i < this->amountOfCities; i++) {
            for (auto j = 0; j < this->amountOfCities; j++) {
                arrayOfMatrixOfCities[i][j] = dist_distancesBetweenCities(gen);
            }
        }
    } else {
        this->amountOfCities = amountOfCities;

        arrayOfMatrixOfCities = new long long int *[this->amountOfCities];
        for (auto i = 0; i < this->amountOfCities; i++) {
            arrayOfMatrixOfCities[i] = new long long int[this->amountOfCities];
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist_distancesBetweenCities(1, maxDistanceBetweenCity);

        for (auto i = 0; i < this->amountOfCities; i++) {
            for (auto j = 0; j < this->amountOfCities; j++) {
                arrayOfMatrixOfCities[i][j] = dist_distancesBetweenCities(gen);
            }
        }
    }
}

void TravellingSalesmanProblem::PrintCitiesForTheTravellingSalesman() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do wyświetlenia.");

    std::cout << "\e[1mProblem\e[0m" << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "Name of TSPLIB file:\t" << fileName << std::endl;
    std::cout << "Graph type:\t\t" << graphType << std::endl;
    std::cout << "Number of cities:\t" << amountOfCities << std::endl;
}

// -------------------------------------------------------------------
// Algorytm zachłanny dla problemu komiwojażera.
// -------------------------------------------------------------------
void TravellingSalesmanProblem::GreedyAlgorithm() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do przeprowadzenia algorytmu problemu komiwojażera.");

    if (optimalWay_Solution != nullptr)
        delete[] optimalWay_Solution;

    setGreedyAlgorithm = true;
    optimalWay_Solution = new int[amountOfCities];

    bool *visitedCities = new bool[amountOfCities];
    for (int i = 0; i < amountOfCities; i++) {
        visitedCities[i] = false;
    }

    length = 0;
    int currentMinLength;

    int nextCity = 0;
    int city = nextCity;
    visitedCities[city] = true;

    optimalWay_Solution[0] = nextCity;

    for (auto j = 0; j < amountOfCities - 1; j++) {
        city = nextCity;
        currentMinLength = INT_MAX;
        for (auto i = 0; i < amountOfCities; i++) {
            if (arrayOfMatrixOfCities[city][i] < currentMinLength && !visitedCities[i]) {
                currentMinLength = arrayOfMatrixOfCities[city][i];
                nextCity = i;
            }
        }
        visitedCities[nextCity] = true;
        optimalWay_Solution[j] = nextCity;
        length += arrayOfMatrixOfCities[city][nextCity];
    }
    optimalWay_Solution[amountOfCities - 1] = 0;
    length += arrayOfMatrixOfCities[optimalWay_Solution[amountOfCities - 2]][0];

    delete[] visitedCities;
}

// -------------------------------------------------------------------
// Algorytm podziału i ograniczeń dla problemu komiwojażera.
// -------------------------------------------------------------------
void TravellingSalesmanProblem::BranchAndBoundAlgorithm() {
    if (arrayOfMatrixOfCities == nullptr)
        throw std::logic_error("Brak miast do przeprowadzenia algorytmu problemu komiwojażera.");

    if (optimalWay_Solution != nullptr)
        delete[] optimalWay_Solution;

    setGreedyAlgorithm = false;
}


void TravellingSalesmanProblem::PrintSolution() {
    std::cout << "\e[1mSolution\e[0m" << std::endl;
    if (setGreedyAlgorithm) {
        std::cout << "\e[1mGreedy Algorithm\e[0m" << std::endl;
    } else {
        std::cout << "\e[1mBranch and Bound Algorithm\e[0m" << std::endl;
    }

    std::cout << "-------------------" << std::endl;
    std::cout << "Length\t= " << length << std::endl;
    std::cout << "Path\t= ";
    if (setGreedyAlgorithm) {
        std::cout << "0 - ";
        for (auto i = 0; i < amountOfCities; i++) {
            if (i == amountOfCities - 1) {
                std::cout << optimalWay_Solution[i] << std::endl;
            } else {
                std::cout << optimalWay_Solution[i] << " - ";
            }
        }
    } else {
        for (auto i = 0; i < amountOfCities; i++) {
            std::cout << optimalWay_Solution[i] << " - ";
        }
        std::cout << "0" << std::endl;
    }
}

