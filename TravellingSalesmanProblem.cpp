//
// Created by mrfarinq on 16.06.17.
//

#include <fstream>
#include <random>
#include <climits>
#include <algorithm>
#include <map>
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
    if (arrayOfMatrixOfCities != nullptr) {
        for (int i = 0; i < amountOfCities; i++)
            delete[] arrayOfMatrixOfCities[i];
        delete[] arrayOfMatrixOfCities;
    }

    amountOfCities = _amountOfCities;

    arrayOfMatrixOfCities = new long long int *[amountOfCities];
    for (int i = 0; i < amountOfCities; i++)
        arrayOfMatrixOfCities[i] = new long long int[amountOfCities];

    for (int i = 0; i < amountOfCities; i++) {
        for (int j = 0; j < amountOfCities; j++) {
            arrayOfMatrixOfCities[i][j] = _cities[i][j];
        }
    }

    fileName = _fileName;
    graphType = _graphType;
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

std::pair<int, long long int> GetMaxValueFromMap(const std::map<int, long long int> &x) {
    using pairtype=std::pair<int, long long int>;
    return *std::max_element(x.begin(), x.end(), [](const pairtype &p1, const pairtype &p2) {
        return p1.second < p2.second;
    });
}

int FindSectionIndex(std::vector<int> _vec, int _amountOfCities, int _indexOfDeletedSection) {
    for (auto i = 0; i < _amountOfCities; i++) {
        if (_vec[i] == _indexOfDeletedSection)
            return i;
    }
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
    optimalWay_Solution = new int[amountOfCities];

    int lowerBound = 0;
    int lowerBoundK1 = 0;
    int lowerBoundK2 = 0;

    std::map<int, int> edgesOfSolution;

    std::map<int, long long int> minimumValueInRow;
    std::map<int, long long int> minimumValueInColumn;

    long long int **activeRoute = new long long int *[amountOfCities];
    long long int **subsetMatrix;
    for (auto i = 0; i < amountOfCities; i++) {
        activeRoute[i] = new long long int[amountOfCities];
    }

    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++) {
            activeRoute[i][j] = arrayOfMatrixOfCities[i][j];
        }
        activeRoute[i][i] = INT_MAX;
    }

    std::vector<int> indexesOfRowsInActiveRoute;
    std::vector<int> indexesOfColumnsInActiveRoute;
    std::vector<int> indexesOfRowsInSubsetK1;
    std::vector<int> indexesOfColumnsInSubsetK1;
    for (auto i = 0; i < amountOfCities; i++) {
        indexesOfRowsInActiveRoute.push_back(i);
        indexesOfColumnsInActiveRoute.push_back(i);
    }

    int amountOfCitiesInActualSubset = amountOfCities;

    //---
    std::cout << "Odległości pomiędzy miastami (macierz wag) oryginalnie: " << std::endl;
    std::cout << "\t";
    for (std::vector<int>::const_iterator i = indexesOfColumnsInActiveRoute.begin();
         i != indexesOfColumnsInActiveRoute.end(); i++) {
        std::cout << *i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
            if (j == 0) {
                if (activeRoute[i][j] < 0) {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t\b" << activeRoute[i][j];
                } else {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t" << activeRoute[i][j];
                }
            } else {
                if (activeRoute[i][j] < 0) {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << activeRoute[i][j];
                } else {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << activeRoute[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie minimalnych wartości dla każdego wiersza.
    minimumValueInRow.clear();
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        minimumValueInRow[indexesOfRowsInActiveRoute[i]] = INT_MAX;
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
            if (activeRoute[i][j] < minimumValueInRow[indexesOfRowsInActiveRoute[i]]) {
                minimumValueInRow[indexesOfRowsInActiveRoute[i]] = activeRoute[i][j];
            }
        }
    }

    //---
    std::cout << "Minimalne wartości w wierszach: " << std::endl;
    for (std::map<int, long long int>::iterator map_iterator = minimumValueInRow.begin();
         map_iterator != minimumValueInRow.end(); ++map_iterator) {
        std::cout << map_iterator->first << " => " << map_iterator->second << std::endl;
    }
    std::cout << std::endl;
    //---

    // Odjęcie minimalnych wartości wiersza od każdego elementu wiersza.
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++)
            if (activeRoute[i][j] != INT_MAX)
                activeRoute[i][j] -= minimumValueInRow[i];
    }

    //---
    std::cout << "Odległości pomiędzy miastami (macierz wag) po odjęciu minimów wierszy: " << std::endl;
    std::cout << "\t";
    for (std::vector<int>::const_iterator i = indexesOfColumnsInActiveRoute.begin();
         i != indexesOfColumnsInActiveRoute.end(); i++) {
        std::cout << *i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
            if (j == 0) {
                if (activeRoute[i][j] < 0) {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t\b" << activeRoute[i][j];
                } else {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t" << activeRoute[i][j];
                }
            } else {
                if (activeRoute[i][j] < 0) {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << activeRoute[i][j];
                } else {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << activeRoute[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie minimalnych wartości dla każdej kolumny.
    minimumValueInColumn.clear();
    for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
        minimumValueInColumn[indexesOfColumnsInActiveRoute[j]] = INT_MAX;
        for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
            if (activeRoute[i][j] < minimumValueInColumn[indexesOfColumnsInActiveRoute[j]]) {
                minimumValueInColumn[indexesOfColumnsInActiveRoute[j]] = activeRoute[i][j];
            }
        }
    }

    //---
    std::cout << std::endl;
    std::cout << "Minimalne wartości w kolumnach: " << std::endl;
    for (std::map<int, long long int>::iterator map_iterator = minimumValueInColumn.begin();
         map_iterator != minimumValueInColumn.end(); ++map_iterator) {
        std::cout << map_iterator->first << " => " << map_iterator->second << std::endl;
    }
    std::cout << std::endl;
    //---

    // Odjęcie minimalnych wartości kolumny od każdego elementu kolumny.
    for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
        for (auto i = 0; i < amountOfCitiesInActualSubset; i++)
            if (activeRoute[i][j] != INT_MAX)
                activeRoute[i][j] -= minimumValueInColumn[j];

    }

    //---
    std::cout << "Odległości pomiędzy miastami (macierz wag) po odjęciu minimów kolumn: " << std::endl;
    std::cout << "\t";
    for (std::vector<int>::const_iterator i = indexesOfRowsInActiveRoute.begin(); i != indexesOfRowsInActiveRoute.end(); i++) {
        std::cout << *i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
            if (j == 0) {
                if (activeRoute[i][j] < 0) {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << indexesOfColumnsInActiveRoute[i] << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfColumnsInActiveRoute[i] << ".\t\b" << activeRoute[i][j];
                } else {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << indexesOfColumnsInActiveRoute[i] << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfColumnsInActiveRoute[i] << ".\t" << activeRoute[i][j];
                }
            } else {
                if (activeRoute[i][j] < 0) {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << activeRoute[i][j];
                } else {
                    if (activeRoute[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << activeRoute[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie parametru lowerBound - suma minimów wierszy + suma minimów kolumn.
    lowerBound = std::accumulate(std::begin(minimumValueInRow), std::end(minimumValueInRow), 0,
                                      [](int value, const std::map<int, int>::value_type &p) {
                                          return value + p.second;
                                      }
    ) + std::accumulate(std::begin(minimumValueInColumn), std::end(minimumValueInColumn), 0,
                        [](int value, const std::map<int, int>::value_type &p) {
                            return value + p.second;
                        }
    );

    //---
    std::cout << std::endl;
    std::cout << "Lower Bound: " << lowerBound << "." << std::endl << std::endl;
    //---

    // Wyznaczenie minimalnych wartości dla każdego wiersza i kolumny (0 uznane za minimum jeżeli wystąpi 2 razy).
    int amountOfZeros;
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        amountOfZeros = 0;
        minimumValueInRow[indexesOfRowsInActiveRoute[i]] = INT_MAX;
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
            if (activeRoute[i][j] == 0)
                amountOfZeros++;
            if (amountOfZeros >= 2) {
                minimumValueInRow[indexesOfRowsInActiveRoute[i]] = 0;
                break;
            }
            if ((activeRoute[i][j] < minimumValueInRow[i]) &&
                (activeRoute[i][j] != 0))
                minimumValueInRow[indexesOfRowsInActiveRoute[i]] = activeRoute[i][j];

        }
    }

    for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
        amountOfZeros = 0;
        minimumValueInColumn[indexesOfColumnsInActiveRoute[j]] = INT_MAX;
        for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
            if (activeRoute[i][j] == 0)
                amountOfZeros++;
            if (amountOfZeros == 2) {
                minimumValueInColumn[indexesOfColumnsInActiveRoute[j]] = 0;
                break;
            }
            if ((activeRoute[i][j] < minimumValueInColumn[j]) &&
                (activeRoute[i][j] != 0))
                minimumValueInColumn[indexesOfColumnsInActiveRoute[j]] = activeRoute[i][j];
        }
    }

    //---
    std::cout << "Minimalne wartości w wierszach (0 uznane za minimum jeżeli wystąpi 2 razy): " << std::endl;
    for (std::map<int, long long int>::iterator map_iterator = minimumValueInRow.begin();
         map_iterator != minimumValueInRow.end(); ++map_iterator) {
        std::cout << map_iterator->first << " => " << map_iterator->second << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Minimalne wartości w kolumnach (0 uznane za minimum jeżeli wystąpi 2 razy): " << std::endl;
    for (std::map<int, long long int>::iterator map_iterator = minimumValueInColumn.begin();
         map_iterator != minimumValueInColumn.end(); ++map_iterator) {
        std::cout << map_iterator->first << " => " << map_iterator->second << std::endl;
    }
    std::cout << std::endl;
    //---

    int maximumValueFromMinimumsOfRows = GetMaxValueFromMap(minimumValueInRow).second;
    int indexOfMaximumValueFromMinimumsOfRows = GetMaxValueFromMap(minimumValueInRow).first;

    int maximumValueFromMinimumsOfColumns = GetMaxValueFromMap(minimumValueInColumn).second;
    int indexOfMaximumValueFromMinimumsOfColumns = GetMaxValueFromMap(minimumValueInColumn).first;

    //---
    if (maximumValueFromMinimumsOfRows > maximumValueFromMinimumsOfColumns)
        std::cout << "Największa wartość minimum - " << maximumValueFromMinimumsOfRows << ", wiersz - "
                  << indexOfMaximumValueFromMinimumsOfRows << "." << std::endl;
    else
        std::cout << "Największa wartość minimum - " << maximumValueFromMinimumsOfColumns << ", kolumna - "
                  << indexOfMaximumValueFromMinimumsOfColumns << "." << std::endl;
    //---

    int indexOfDeletedRow;
    int indexOfDeletedColumn;
    int indexOfDeletedSection;

    long long int **K2 = new long long int *[amountOfCitiesInActualSubset];
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        K2[i] = new long long int[amountOfCitiesInActualSubset];
    }

    for (int i = 0; i < amountOfCitiesInActualSubset; i++) {
        for (int j = 0; j < amountOfCitiesInActualSubset; j++) {
            K2[i][j] = activeRoute[i][j];
        }
    }

    if (maximumValueFromMinimumsOfRows > maximumValueFromMinimumsOfColumns) {
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
            if (activeRoute[indexOfMaximumValueFromMinimumsOfRows][j] == 0) {
                indexOfDeletedSection = j;
                indexOfDeletedRow = indexOfMaximumValueFromMinimumsOfRows;
                indexOfDeletedColumn = indexOfDeletedSection;
                lowerBoundK2 = lowerBound + maximumValueFromMinimumsOfRows;
                break;
            }
        }
        K2[indexOfMaximumValueFromMinimumsOfRows][indexOfDeletedSection]=INT_MAX;
        activeRoute[indexOfDeletedSection][indexOfMaximumValueFromMinimumsOfRows] = INT_MAX;
    } else {

        for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
            if (activeRoute[i][indexOfMaximumValueFromMinimumsOfColumns] == 0) {
                indexOfDeletedSection = i;
                indexOfDeletedRow = indexOfDeletedSection;
                indexOfDeletedColumn = indexOfMaximumValueFromMinimumsOfColumns;
                lowerBoundK2 = lowerBound + maximumValueFromMinimumsOfColumns;
                break;
            }
        }
        K2[indexOfDeletedSection][indexOfMaximumValueFromMinimumsOfColumns]=INT_MAX;
        activeRoute[indexOfMaximumValueFromMinimumsOfColumns][indexOfDeletedSection] = INT_MAX;
    }


    //Usunięcie indeksów usuwanych wierszy i kolumn.
    indexesOfRowsInSubsetK1.clear();
    indexesOfColumnsInSubsetK1.clear();
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        indexesOfRowsInSubsetK1.push_back(indexesOfRowsInActiveRoute[i]);
        indexesOfColumnsInSubsetK1.push_back(indexesOfColumnsInActiveRoute[i]);
    }
    indexesOfRowsInSubsetK1.erase(indexesOfRowsInSubsetK1.begin() +
                                FindSectionIndex(indexesOfRowsInSubsetK1, amountOfCitiesInActualSubset,
                                                 indexOfDeletedRow));
    indexesOfColumnsInSubsetK1.erase(indexesOfColumnsInSubsetK1.begin() +
                                   FindSectionIndex(indexesOfColumnsInSubsetK1, amountOfCitiesInActualSubset,
                                                    indexOfDeletedColumn));

    //---
    std::cout << std::endl;
    std::cout << "Wypisanie aktualnych wierszy K1: " << std::endl;
    for (std::vector<int>::const_iterator i = indexesOfRowsInSubsetK1.begin(); i != indexesOfRowsInSubsetK1.end(); i++) {
        std::cout << *i << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Wypisanie aktualnych kolumn K1: " << std::endl;
    for (std::vector<int>::const_iterator i = indexesOfColumnsInSubsetK1.begin();
         i != indexesOfColumnsInSubsetK1.end(); i++) {
        std::cout << *i << std::endl;
    }
    std::cout << std::endl;
    //---

    long long int **K1 = new long long int *[amountOfCitiesInActualSubset - 1];
    for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++) {
        K1[i] = new long long int[amountOfCitiesInActualSubset - 1];
    }

    //Redukcja macierzy K1 - usunięcie wiersza.
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++)
        if (i != indexOfDeletedRow) {
            for (auto j = 0; j < amountOfCitiesInActualSubset; j++)
                if (j != indexOfDeletedColumn) {
                    if (i < indexOfDeletedRow && j < indexOfDeletedColumn)
                        K1[i][j] = activeRoute[i][j];
                    else if (i < indexOfDeletedRow && j > indexOfDeletedColumn)
                        K1[i][j - 1] = activeRoute[i][j];
                    else if (i > indexOfDeletedRow && j < indexOfDeletedColumn)
                        K1[i - 1][j] = activeRoute[i][j];
                    else
                        K1[i - 1][j - 1] = activeRoute[i][j];
                }
        }

    //---
    std::cout << "Odległości pomiędzy miastami (zredukowana macierz wag K1) z dodaną blokadą: " << std::endl;
    std::cout << "\t";
    for (std::vector<int>::const_iterator i = indexesOfColumnsInSubsetK1.begin();
         i != indexesOfColumnsInSubsetK1.end(); i++) {
        std::cout << *i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset - 1; j++) {
            if (j == 0) {
                if (K1[i][j] < 0) {
                    if (K1[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t\b" << K1[i][j];
                } else {
                    if (K1[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t" << K1[i][j];
                }
            } else {
                if (K1[i][j] < 0) {
                    if (K1[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << K1[i][j];
                } else {
                    if (K1[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << K1[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie minimalnych wartości dla każdego wiersza K1.
    minimumValueInRow.clear();
    for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++) {
        minimumValueInRow[indexesOfRowsInSubsetK1[i]] = INT_MAX;
        for (auto j = 0; j < amountOfCitiesInActualSubset - 1; j++) {
            if (K1[i][j] < minimumValueInRow[indexesOfRowsInSubsetK1[i]]) {
                minimumValueInRow[indexesOfRowsInSubsetK1[i]] = K1[i][j];
            }
        }
    }

    //---
    std::cout << "Minimalne wartości w wierszach K1: " << std::endl;
    for (std::map<int, long long int>::iterator map_iterator = minimumValueInRow.begin();
         map_iterator != minimumValueInRow.end(); ++map_iterator) {
        std::cout << map_iterator->first << " => " << map_iterator->second << std::endl;
    }
    std::cout << std::endl;
    //---

    // Odjęcie minimalnych wartości wiersza od każdego elementu wiersza.
    for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset - 1; j++)
            if (K1[i][j] != INT_MAX)
                K1[i][j] -= minimumValueInRow[i];
    }

    //---
    std::cout
            << "Odległości pomiędzy miastami (zredukowana macierz wag K1) z dodaną blokadą po 1. etapie standaryzacji: "
            << std::endl;
    std::cout << "\t";
    for (std::vector<int>::const_iterator i = indexesOfColumnsInSubsetK1.begin();
         i != indexesOfColumnsInSubsetK1.end(); i++) {
        std::cout << *i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset - 1; j++) {
            if (j == 0) {
                if (K1[i][j] < 0) {
                    if (K1[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t\b" << K1[i][j];
                } else {
                    if (K1[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t" << K1[i][j];
                }
            } else {
                if (K1[i][j] < 0) {
                    if (K1[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << K1[i][j];
                } else {
                    if (K1[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << K1[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie minimalnych wartości dla każdej kolumny K1.
    minimumValueInColumn.clear();
    for (auto j = 0; j < amountOfCitiesInActualSubset - 1; j++) {
        minimumValueInColumn[indexesOfColumnsInSubsetK1[j]] = INT_MAX;
        for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++) {
            if (K1[i][j] < minimumValueInColumn[indexesOfColumnsInSubsetK1[j]]) {
                minimumValueInColumn[indexesOfColumnsInSubsetK1[j]] = K1[i][j];
            }
        }
    }

    //---
    std::cout << std::endl;
    std::cout << "Minimalne wartości w kolumnach K1: " << std::endl;
    for (std::map<int, long long int>::iterator map_iterator = minimumValueInColumn.begin();
         map_iterator != minimumValueInColumn.end(); ++map_iterator) {
        std::cout << map_iterator->first << " => " << map_iterator->second << std::endl;
    }
    std::cout << std::endl;
    //---

    // Odjęcie minimalnych wartości kolumny od każdego elementu kolumny.
    for (auto j = 0; j < amountOfCitiesInActualSubset - 1; j++) {
        for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++)
            if (K1[i][j] != INT_MAX)
                K1[i][j] -= minimumValueInColumn[j];

    }

    //---
    std::cout
            << "Odległości pomiędzy miastami (zredukowana macierz wag K1) z dodaną blokadą po 2. etapie standaryzacji: "
            << std::endl;
    std::cout << "\t";
    for (std::vector<int>::const_iterator i = indexesOfColumnsInSubsetK1.begin();
         i != indexesOfColumnsInSubsetK1.end(); i++) {
        std::cout << *i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCitiesInActualSubset - 1; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset - 1; j++) {
            if (j == 0) {
                if (K1[i][j] < 0) {
                    if (K1[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t\b" << K1[i][j];
                } else {
                    if (K1[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInSubsetK1[i] << ".\t" << K1[i][j];
                }
            } else {
                if (K1[i][j] < 0) {
                    if (K1[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << K1[i][j];
                } else {
                    if (K1[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << K1[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie parametru lowerBound dla K1 - poprzedni lowerBound + suma minimów wierszy + suma minimów kolumn.
    lowerBoundK1 = lowerBound + std::accumulate(std::begin(minimumValueInRow), std::end(minimumValueInRow), 0,
                                                [](int value, const std::map<int, int>::value_type &p) {
                                                    return value + p.second;
                                                }
    ) + std::accumulate(std::begin(minimumValueInColumn), std::end(minimumValueInColumn), 0,
                        [](int value, const std::map<int, int>::value_type &p) {
                            return value + p.second;
                        }
    );

    //---
    std::cout << "Lower Bound K1: " << lowerBoundK1 << ", Lower Bound K2: " << lowerBoundK2 << "." << std::endl
              << std::endl;
    //---

    //---
    std::cout
            << "Odległości pomiędzy miastami (macierz wag K2) z dodaną blokadą: "
            << std::endl;
    std::cout << "\t";
    for (std::vector<int>::const_iterator i = indexesOfColumnsInActiveRoute.begin();
         i != indexesOfColumnsInActiveRoute.end(); i++) {
        std::cout << *i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCitiesInActualSubset; i++) {
        for (auto j = 0; j < amountOfCitiesInActualSubset; j++) {
            if (j == 0) {
                if (K2[i][j] < 0) {
                    if (K2[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t\b" << K2[i][j];
                } else {
                    if (K2[i][j] == INT_MAX)
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << indexesOfRowsInActiveRoute[i] << ".\t" << K2[i][j];
                }
            } else {
                if (K2[i][j] < 0) {
                    if (K2[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << K2[i][j];
                } else {
                    if (K2[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << K2[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    if (lowerBoundK1>lowerBoundK2)
    {

    }
    else
    {

    }

    //edgesOfSolution.insert(indexOfDeletedRow, indexOfDeletedColumn);
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

