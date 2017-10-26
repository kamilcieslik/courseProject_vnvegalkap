//
// Created by mrfarinq on 16.06.17.
//

#ifndef SDIZO_3_TRAVELLINGSALESMANPROBLEM_H
#define SDIZO_3_TRAVELLINGSALESMANPROBLEM_H

#include <iostream>

class TravellingSalesmanProblem {
private:
    int amountOfCities;
    long long int **arrayOfMatrixOfCities;
    int *optimalWay_Solution;
    int length;
    bool setGreedyAlgorithm;
    std::string fileName;
    std::string graphType;

public:
    TravellingSalesmanProblem();

    ~TravellingSalesmanProblem();

    void DeleteTravellingSalesman();

    void LoadArrayOfMatrixOfCities(long long int **_cities, int _amountOfCities, std::string _fileName,
                                   std::string _graphType);

    void GenerateRandomCities(int amountOfCities = 0, int maxDistanceBetweenCity = 99);

    void PrintCitiesForTheTravellingSalesman();

    void GreedyAlgorithm();

    void BranchAndBoundAlgorithm();

    void PrintSolution();
};


#endif //SDIZO_3_TRAVELLINGSALESMANPROBLEM_H
