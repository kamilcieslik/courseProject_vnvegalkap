//
// Created by mrfarinq on 26.10.17.
//

#ifndef COURSEPROJECT_SDZ3_TSPLIB_PARSER_H
#define COURSEPROJECT_SDZ3_TSPLIB_PARSER_H


#include <vector>
#include <string>

class TSPLIB_Parser {
private:
    const char delimiter = ':';
    std::string name;
    std::string type;
    std::string comment;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::vector<long long int> numbers;
    std::vector<long long int> optimalTour;
    int cost;

public:
    long long int **cities = nullptr;// Interpreted Matrix, usable for the Little algorithm
    int dimension;

    bool checkParameter(std::string, std::string);

    std::string trim(std::string);

    bool readProblem(std::ifstream &);

    void printSolution();

    void writeSolution(std::ofstream &);

    bool fillMatrix();

    void euclidesMatrix();

    void pseudoEuclidesMatrix();

    void fullMatrix();

    void upperRow();

    void lowerRow();

    void upperDiagRow();

    void lowerDiagRow();

    TSPLIB_Parser(std::ifstream &);

    ~TSPLIB_Parser();
};


#endif //COURSEPROJECT_SDZ3_TSPLIB_PARSER_H
