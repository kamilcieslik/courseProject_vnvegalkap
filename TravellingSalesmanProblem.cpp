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

    std::map<int, int> edgesOfSolution;

    std::map<int, long long int> minimumValueInRow;
    std::map<int, long long int> minimumValueInColumn;

    long long int **arrayOfBranchAndBoundMatrixOfCities = new long long int *[amountOfCities];
    for (auto i = 0; i < amountOfCities; i++)
        arrayOfBranchAndBoundMatrixOfCities[i] = new long long int[amountOfCities];

    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++)
            arrayOfBranchAndBoundMatrixOfCities[i][j] = arrayOfMatrixOfCities[i][j];
        arrayOfBranchAndBoundMatrixOfCities[i][i] = INT_MAX;
    }

    int amountOfCitiesInBranchAndBoundMatrix = amountOfCities;

    //---
    std::cout << "Odległości pomiędzy miastami (macierz wag) oryginalnie: " << std::endl;
    std::cout << "\t";
    for (auto i = 0; i < amountOfCities; i++) {
        std::cout << i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++) {
            if (j == 0) {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            } else {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie minimalnych wartości dla każdego wiersza.
    for (auto i = 0; i < amountOfCitiesInBranchAndBoundMatrix; i++) {
        minimumValueInRow[i] = INT_MAX;
        for (auto j = 0; j < amountOfCitiesInBranchAndBoundMatrix; j++) {
            if (arrayOfBranchAndBoundMatrixOfCities[i][j] < minimumValueInRow[i]) {
                minimumValueInRow[i] = arrayOfBranchAndBoundMatrixOfCities[i][j];
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
    for (auto i = 0; i < amountOfCitiesInBranchAndBoundMatrix; i++) {
        for (auto j = 0; j < amountOfCitiesInBranchAndBoundMatrix; j++)
            if (arrayOfBranchAndBoundMatrixOfCities[i][j] != INT_MAX)
                arrayOfBranchAndBoundMatrixOfCities[i][j] -= minimumValueInRow[i];
    }

    //---
    std::cout << "Odległości pomiędzy miastami (macierz wag) po odjęciu minimów wierszy: " << std::endl;
    std::cout << "\t";
    for (auto i = 0; i < amountOfCities; i++) {
        std::cout << i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++) {
            if (j == 0) {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            } else {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie minimalnych wartości dla każdej kolumny.
    for (auto j = 0; j < amountOfCitiesInBranchAndBoundMatrix; j++) {
        minimumValueInColumn[j] = INT_MAX;
        for (auto i = 0; i < amountOfCitiesInBranchAndBoundMatrix; i++) {
            if (arrayOfBranchAndBoundMatrixOfCities[i][j] < minimumValueInColumn[j]) {
                minimumValueInColumn[j] = arrayOfBranchAndBoundMatrixOfCities[i][j];
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
    for (auto j = 0; j < amountOfCitiesInBranchAndBoundMatrix; j++) {
        for (auto i = 0; i < amountOfCitiesInBranchAndBoundMatrix; i++)
            if (arrayOfBranchAndBoundMatrixOfCities[i][j] != INT_MAX)
                arrayOfBranchAndBoundMatrixOfCities[i][j] -= minimumValueInColumn[j];

    }

    //---
    std::cout << "Odległości pomiędzy miastami (macierz wag) po odjęciu minimów kolumn: " << std::endl;
    std::cout << "\t";
    for (auto i = 0; i < amountOfCities; i++) {
        std::cout << i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++) {
            if (j == 0) {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            } else {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    // Wyznaczenie parametru lowerBound - suma minimów wierszy + suma minimów kolumn.
    auto lowerBound = std::accumulate(std::begin(minimumValueInRow), std::end(minimumValueInRow), 0,
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
    for (auto i = 0; i < amountOfCitiesInBranchAndBoundMatrix; i++) {
        amountOfZeros = 0;
        minimumValueInRow[i] = INT_MAX;
        for (auto j = 0; j < amountOfCitiesInBranchAndBoundMatrix; j++) {
            if (arrayOfBranchAndBoundMatrixOfCities[i][j] == 0)
                amountOfZeros++;
            if (amountOfZeros >= 2) {
                minimumValueInRow[i] = 0;
                break;
            }
            if ((arrayOfBranchAndBoundMatrixOfCities[i][j] < minimumValueInRow[i]) &&
                (arrayOfBranchAndBoundMatrixOfCities[i][j] != 0))
                minimumValueInRow[i] = arrayOfBranchAndBoundMatrixOfCities[i][j];

        }
    }

    for (auto j = 0; j < amountOfCitiesInBranchAndBoundMatrix; j++) {
        amountOfZeros = 0;
        minimumValueInColumn[j] = INT_MAX;
        for (auto i = 0; i < amountOfCitiesInBranchAndBoundMatrix; i++) {
            if (arrayOfBranchAndBoundMatrixOfCities[i][j] == 0)
                amountOfZeros++;
            if (amountOfZeros == 2) {
                minimumValueInColumn[j] = 0;
                break;
            }
            if ((arrayOfBranchAndBoundMatrixOfCities[i][j] < minimumValueInColumn[j]) &&
                (arrayOfBranchAndBoundMatrixOfCities[i][j] != 0))
                minimumValueInColumn[j] = arrayOfBranchAndBoundMatrixOfCities[i][j];
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

    if (maximumValueFromMinimumsOfRows > maximumValueFromMinimumsOfColumns) {
        for (auto j = 0; j < amountOfCitiesInBranchAndBoundMatrix; j++) {
            if (arrayOfBranchAndBoundMatrixOfCities[indexOfMaximumValueFromMinimumsOfRows][j] == 0) {
                indexOfDeletedSection = j;
                indexOfDeletedRow=indexOfMaximumValueFromMinimumsOfRows;
                indexOfDeletedColumn=indexOfDeletedSection;
                break;
            }
        }
        arrayOfBranchAndBoundMatrixOfCities[indexOfDeletedSection][indexOfMaximumValueFromMinimumsOfRows] = INT_MAX;
    } else {

        for (auto i = 0; i < amountOfCitiesInBranchAndBoundMatrix; i++) {
            if (arrayOfBranchAndBoundMatrixOfCities[i][indexOfMaximumValueFromMinimumsOfColumns] == 0) {
                indexOfDeletedSection = i;
                indexOfDeletedRow=indexOfDeletedSection;
                indexOfDeletedColumn=indexOfMaximumValueFromMinimumsOfRows;
                break;
            }
        }
        arrayOfBranchAndBoundMatrixOfCities[indexOfMaximumValueFromMinimumsOfColumns][indexOfDeletedSection] = INT_MAX;
    }

    //---
    std::cout << "Odległości pomiędzy miastami (macierz wag) z dodaną blokadą: " << std::endl;
    std::cout << "\t";
    for (auto i = 0; i < amountOfCities; i++) {
        std::cout << i << ".\t";
    }
    std::cout << "\v" << std::endl;
    for (auto i = 0; i < amountOfCities; i++) {
        for (auto j = 0; j < amountOfCities; j++) {
            if (j == 0) {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << i << ".\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << i << ".\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            } else {
                if (arrayOfBranchAndBoundMatrixOfCities[i][j] < 0) {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t\b" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t\b" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                } else {
                    if (arrayOfBranchAndBoundMatrixOfCities[i][j] == INT_MAX)
                        std::cout << "\t" << "\e[1mINF\e[0m";
                    else
                        std::cout << "\t" << arrayOfBranchAndBoundMatrixOfCities[i][j];
                }
            }
        }
        std::cout << "\v" << std::endl;
    }
    //---

    edgesOfSolution.insert(indexOfDeletedRow, indexOfDeletedColumn);
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

