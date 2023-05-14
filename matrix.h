#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

template<typename T>
class Matrix {
    int rows, cols;
    std::vector<std::vector<T>> matrix = {};

public:
    Matrix();
    Matrix(int, int);
    explicit Matrix(const std::string&);
    Matrix(const Matrix<T>&);
    Matrix(std::vector<std::vector<T>>);
   ~Matrix() = default;


    T         det() const;
    Matrix<T> sub_matrix(int, int) const;
    void      transpose();

    void set_matrix_value(int i, int j, T a) 
        { matrix[i][j] = a; }

    std::vector<std::vector<T>>  get_matrix() const
        { return this->matrix; }
        
    int get_rows() const { return rows; }
    int get_cols() const { return cols; }

    Matrix<T> operator*(const Matrix<T>&) const;
    Matrix<T> operator*(const T&) const;

    Matrix<T> operator+(const Matrix<T>&) const;
    Matrix<T> operator-(const Matrix<T>&) const;

    Matrix<float> operator!() const;

    bool operator==(const Matrix<T>&) const;
    bool operator==(const int&) const;

    static Matrix<int> eye(int, int);
    static Matrix<int> zero(int, int);
};
///////////////////////////////////////////////////////////////////////////
//#######################################################################//
template<typename T>
Matrix<T>::Matrix() {
    std::cin >> this->rows >> this->cols;
    for (int i = 0; i < this->rows; ++i) {
        std::vector<T> row = {};
        for (int j = 0; j < this->cols; ++j) {
            T value;
            std::cin >> value;
            row.push_back(value);
        }
        this->matrix.push_back(row);
    }
}

template<typename T>
Matrix<T>::Matrix(const std::string& filename) {
    std::ifstream file;
    file.open(filename);

    file >> this->rows >> this->cols;
    for (int i = 0; i < this->rows; ++i) {
        std::vector<T> row = {};
        for (int j = 0; j < this->cols; ++j) {
            T value;
            file >> value;
            row.push_back(value);
        }
        this->matrix.push_back(row);
    }

    file.close();
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& m)
    : rows(m.rows), cols(m.cols), matrix(m.matrix) {}

template<typename T>
Matrix<T>::Matrix(std::vector<std::vector<T>> m)
    : rows(m.size()), cols(m[0].size()), matrix(m) {}

template<typename T>
Matrix<T>::Matrix(int _rows, int _cols) {
    rows = _rows; cols = _cols;
    for (int i = 0; i < rows; ++i) {
        std::vector<T> row = {};
        for (int j = 0; j < cols; ++j) {
            row.push_back(0);
        }
        matrix.push_back(row);
    }
}
//#####################################################//

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& m2) const {
    std::vector<std::vector<T>> result;
    for (int i = 0; i < this->rows; ++i) {
        std::vector<T> row = {};
        for (int j = 0; j < m2.cols; ++j) {
            T elem = 0;
            for (int k = 0; k < m2.rows; ++k)
                elem += this->matrix[i][k] * m2.matrix[k][j];
            row.push_back(elem);
        }
        result.push_back(row);
    }
    Matrix<T> ret = result;
    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& a) const {
    Matrix<T> ret = this->matrix;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            ret.matrix[i][j] *= a;
    return ret;
}

template<typename T>
Matrix<T> operator*(const T& a, const Matrix<T>& m) {
    Matrix<T> ret = m * a;
    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& m) const {
    Matrix<T> ret = matrix;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            ret.matrix[i][j] += m.matrix[i][j];
    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& m) const {
    Matrix<T> ret = matrix;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            ret.matrix[i][j] -= m.matrix[i][j];
    return ret;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& m) const {
    if (rows != m.rows || cols != m.cols)
        return false;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (matrix[i][j] != m.matrix[i][j])
                return false;
    return true;
}

template<typename T>
bool Matrix<T>::operator==(const int& a) const {
    if (a != 1 && a != 0)
        exit(1);
    if (a == 1)
            return ((*this) == Matrix<int>::eye(this->rows, this->cols));
    if (a == 0)
        return ((*this) == Matrix<int>::zero(this->rows, this->cols));
    return false;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& m) {
    for (int i = 0; i < m.get_rows(); ++i) {
        stream << ">  ";
        for (int j = 0; j < m.get_cols(); ++j)
            stream << m.get_matrix()[i][j] << "\t";
        stream << "\n";
    }
    return stream;
}
//########//######//########//##############//#########//############//

template<typename T>
Matrix<T> Matrix<T>::sub_matrix(int _row,int _col) const {
    std::vector<std::vector<T>> res = {};
    for (int i = 0; i < rows; ++i) {
        if (i == _row)
            continue;
        std::vector<T> row = {};
        for (int j = 0; j < cols; ++j) {
            if (j == _col)
                continue;
            row.push_back(matrix[i][j]);
        }
        res.push_back(row);
    }
    Matrix<T> ret = res;
    return ret;
}


template<typename T>
T Matrix<T>::det() const {
    if (rows != cols)
        exit(1);
    if (rows == 1)
        return matrix[0][0];
    else if (rows == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    T ret = 0;
    int sign = 1;
    for (int i = 0; i < cols; ++i) {
        Matrix<T> sub = this->sub_matrix(0, i);
        ret += sign * matrix[0][i] * sub.det();
        sign = -sign;
    }
    return ret;
}

template<typename T>
void Matrix<T>::transpose() {
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < i; ++j) {
            T tmp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = tmp;
        }
}

template<typename T>
Matrix<float> Matrix<T>::operator!() const {
    T _det = this->det();
    Matrix<float> alg_adj = this->matrix;
    int sign = 1;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Matrix<T> sub = this->sub_matrix(i, j);
            if ((i + j) % 2 == 0)
                sign = 1;
            else 
                sign = -1;
            alg_adj.set_matrix_value(i, j, sign * sub.det());
        }
    }

    alg_adj.transpose();
    alg_adj = alg_adj * (1/_det);
    return alg_adj;
}
///////////////////////////////////////////////////////////

template<typename T>
Matrix<int> Matrix<T>::eye(int _rows, int _cols) {
    Matrix<int> ret(_rows, _cols);
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            if (i == j)
                ret.set_matrix_value(i, j, 1);
        }
    }
    return ret;
}

template<typename T>
Matrix<int> Matrix<T>::zero(int _rows, int _cols) {
    Matrix<int> ret(_rows, _cols);
    return ret;
}

#endif
