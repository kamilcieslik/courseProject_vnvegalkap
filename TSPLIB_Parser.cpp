//
// Created by mrfarinq on 26.10.17.
//

#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "TSPLIB_Parser.h"


TSPLIB_Parser::TSPLIB_Parser(std::ifstream &inputFile) {
    if (readProblem(inputFile)) {
        // printSolution();
    }
}

TSPLIB_Parser::~TSPLIB_Parser() {
    if(cities!=nullptr){
        for (int i = 0; i < dimension; i++)
            delete[] cities[i];
        delete[] cities;
    }
    cities=nullptr;
}

bool TSPLIB_Parser::readProblem(std::ifstream &inputFile) {
    std::string line;
    bool isMatrixType = false;
    bool isCoordinatesType = false;
    while (inputFile) {
        getline(inputFile, line);
        if (line == "EOF" || line == "DISPLAY_DATA_SECTION") {
            break;
        }
        if (line.find(delimiter) != std::string::npos) {
            std::string parameter = line.substr(0, line.find(delimiter));
            std::string value = line.substr(line.find(delimiter) + 1, line.npos);

            if (!checkParameter(trim(parameter), trim(value))) {
                return false;
            }
        }
        if (isMatrixType) {
            std::stringstream stream(line);
            int n;
            while (stream >> n) {
                this->numbers.push_back((long long int &&) n);
            }
        }

        if (isCoordinatesType) {
            std::stringstream stream(line);
            int n;
            stream >> n;
            while (stream >> n) {
                this->numbers.push_back((long long int &&) n);
            }

        }
        if (line == "EDGE_WEIGHT_SECTION") {
            isMatrixType = true;
        }
        if (line == "NODE_COORD_SECTION") {
            isCoordinatesType = true;
        }
    }

    fillMatrix();


    return true;
}

bool TSPLIB_Parser::checkParameter(std::string keyword, std::string value) {
    if (keyword == "NAME") {
        this->name = value;
    } else if (keyword == "TYPE") {
        if ((value == "TSP") || (value == "ATSP")) {
            this->type = value;
        } else {
            std::cout << "Parametr - " << keyword << " niewspierany." << std::endl;
            return 0;
        }
    } else if (keyword == "COMMENT")
        this->comment = value;
    else if (keyword == "DIMENSION")
        this->dimension = stoi(value);
    else if (keyword == "EDGE_WEIGHT_TYPE") {
        if (value == "EXPLICIT")
            this->edgeWeightType = value;
        else if (value == "EUC_2D")
            this->edgeWeightType = value;
        else if (value == "ATT")
            this->edgeWeightFormat = value;
        else {
            std::cout << "Parametr - " << keyword << " niewspierany." << std::endl;
            return 0;
        }
    } else if (keyword == "EDGE_WEIGHT_FORMAT") {
        if ((value == "FULL_MATRIX") ||
            (value == "UPPER_ROW") ||
            (value == "LOWER_ROW") ||
            (value == "UPPER_DIAG_ROW") ||
            (value == "LOWER_DIAG_ROW") ||
            (value == "UPPER_COL") ||
            (value == "LOWER_COL") ||
            (value == "UPPER_DIAG_COL") ||
            (value == "LOWER_DIAG_COL"))
            this->edgeWeightFormat = value;
        else {
            std::cout << "Parametr - " << keyword << " niewspierany." << std::endl;
            return 0;
        }
    } else if (keyword == "DISPLAY_DATA_TYPE") {

    } else {
        std::cout << "Nieznany parametr - " << keyword << "." << std::endl;
        return false;
    }
    return true;
}

std::string TSPLIB_Parser::trim(std::string s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

void TSPLIB_Parser::printSolution() {
    time_t now = time(0);
    tm *localtm = localtime(&now);

    std::cout << "NAME : " << this->name << "." << this->optimalTour.size() << ".tour" << std::endl;
    std::cout << "COMMENT : Lenght = " << this->cost << ". Found by John D.C. Little " << asctime(localtm);
    std::cout << "TYPE : TOUR" << std::endl;
    std::cout << "DIMENSION : " << this->optimalTour.size() << std::endl;
    std::cout << "TOUR_SECTION" << std::endl;

    for (long long int & i : this->optimalTour) {
        std::cout << i << " ";
    }
    std::cout << "-1" << std::endl;
}

void TSPLIB_Parser::writeSolution(std::ofstream &outputFile) {
    time_t now = time(0);
    tm *localtm = localtime(&now);

    outputFile << "NAME : " << this->name << "." << this->optimalTour.size() << ".tour" << std::endl;
    outputFile << "COMMENT : Lenght = " << this->cost << ". Found by John D.C. Little " << asctime(localtm);
    outputFile << "TYPE : TOUR" << std::endl;
    outputFile << "DIMENSION : " << this->optimalTour.size() << std::endl;
    outputFile << "TOUR_SECTION" << std::endl;

    for (long long int & i : this->optimalTour) {
        outputFile << i << std::endl;
    }
    outputFile << "-1" << std::endl;
    outputFile << "EOF";
}

bool TSPLIB_Parser::fillMatrix() {
    if (cities != nullptr) {
        for (int i = 0; i < dimension; i++)
            delete[] cities[i];
        delete[] cities;
    }
    cities = new long long int *[dimension];
    for (int i = 0; i < dimension; i++)
        cities[i] = new long long int[dimension];

    if (edgeWeightType == "EUC_2D")
        euclidesMatrix();

    if (edgeWeightFormat=="ATT")
        pseudoEuclidesMatrix();

    if (edgeWeightFormat == "FULL_MATRIX")
        fullMatrix();

    if ((edgeWeightFormat == "UPPER_ROW") || (edgeWeightFormat == "LOWER_COL"))
        upperRow();

    if ((edgeWeightFormat == "LOWER_ROW") || (edgeWeightFormat == "UPPER_COL"))
        lowerRow();

    if ((edgeWeightFormat == "UPPER_DIAG_ROW") || (edgeWeightFormat == "LOWER_DIAG_COL"))
        upperDiagRow();

    if ((edgeWeightFormat == "LOWER_DIAG_ROW") || (edgeWeightFormat == "UPPER_DIAG_COL"))
        lowerDiagRow();


    return true;
}

void TSPLIB_Parser::fullMatrix() {
    for (int i = 0; i < this->dimension; i++) {
        for (int j = 0; j < this->dimension; j++) {
            if (i != j) {
                cities[i][j] = this->numbers[i * this->dimension + j];
            }
        }
    }
}

void TSPLIB_Parser::upperRow() {
    int counter = 0;
    for (int i = 0; i < this->dimension - 1; i++) {
        for (int j = i + 1; j < this->dimension; j++) {
            cities[i][j] = this->numbers[counter];
            cities[j][i] = this->numbers[counter];
            counter++;
        }
    }
}

void TSPLIB_Parser::lowerRow() {
    int counter = 0;
    for (int i = 1; i < this->dimension; i++) {
        for (int j = 0; j < i; j++) {
            cities[i][j] = this->numbers[counter];
            cities[j][i] = this->numbers[counter];
            counter++;
        }
    }
}

void TSPLIB_Parser::upperDiagRow() {
    int counter = 0;
    for (int i = 0; i < this->dimension; i++) {
        for (int j = i; j < this->dimension; j++) {
            if (i != j) {
                cities[i][j] = this->numbers[counter];
                cities[j][i] = this->numbers[counter];
            }
            counter++;
        }
    }
}

void TSPLIB_Parser::lowerDiagRow() {
    int counter = 0;
    for (int i = 0; i < this->dimension; i++) {
        for (int j = 0; j < i + 1; j++) {
            if (j != i) {
                cities[i][j] = this->numbers[counter];
                cities[j][i] = this->numbers[counter];
            }
            counter++;
        }
    }
}

void TSPLIB_Parser::euclidesMatrix() {
    int counter=0;
    int counter2=0;
    for (int i = 0; i < (2*dimension); i=i+2) {
        for (int j = 0; j <(2*dimension) ; j=j+2) {

            cities[counter][counter2] =  (long long int)round(sqrt( (numbers[i]-numbers[j])*(numbers[i]-numbers[j]) + (numbers[i+1]-numbers[j+1])*(numbers[i+1]-numbers[j+1]) ));
            counter2++;
        }
        std::cout<<"kurwa"<<std::endl;
        counter2=0;
        counter ++;
    }
    std::cout<<std::endl;
    std::cout<<counter;
    std::cout<<"wyszlo git";
}

void TSPLIB_Parser::pseudoEuclidesMatrix()
{
    double rij;
    long long int tij;
    int counter=0;
    int counter2=0;
    for (int i = 0; i < (2*dimension); i=i+2) {
        for (int j = 0; j <(2*dimension) ; j=j+2) {

            rij =  sqrt( ((numbers[i]-numbers[j])*(numbers[i]-numbers[j]) + (numbers[i+1]-numbers[j+1])*(numbers[i+1]-numbers[j+1]))/10.0 );
            tij = (long long int)round(rij);
            if (tij<rij) cities[counter][counter2] = tij + 1;
            else cities[counter][counter2] = tij;

            counter2++;
        }
        counter2=0;
        counter ++;
    }
    std::cout<<std::endl;
    std::cout<<counter;
}
