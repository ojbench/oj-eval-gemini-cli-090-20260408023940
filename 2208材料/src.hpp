#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"
#include <vector>
#include <utility>

class matrix {
private:

    // m行n列的矩阵，用动态二维数组存储，每个元素是分数类实例
    int m, n;
    fraction **data;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************
    friend class resistive_network;

public:

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

    // 默认构造函数
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // TODO: 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
    matrix(int m_, int n_) {
        if (m_ <= 0 || n_ <= 0) {
            m = n = 0;
            data = nullptr;
        } else {
            m = m_;
            n = n_;
            data = new fraction*[m];
            for (int i = 0; i < m; ++i) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; ++j) {
                    data[i][j] = fraction(0);
                }
            }
        }
    }

    // TODO: 拷贝构造函数，构建与 obj 完全相同的矩阵。
    matrix(const matrix &obj) {
        m = obj.m;
        n = obj.n;
        if (m == 0 || n == 0) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; ++i) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; ++j) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
    }

    // TODO: 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept {
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    // TODO: 析构函数。
    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    matrix &operator=(matrix &&obj) noexcept {
        if (this == &obj) return *this;
        if (data != nullptr) {
            for (int i = 0; i < m; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
        return *this;
    }

    // TODO: 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        if (data != nullptr) {
            for (int i = 0; i < m; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }
        m = obj.m;
        n = obj.n;
        if (m == 0 || n == 0) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; ++i) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; ++j) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
        return *this;
    }

    // TODO: 重载括号，返回矩阵的第i行(1-based)、第j列(0-based)的元素的引用。如果 i、j 不合法，抛出 matrix_error 错误。
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n) {
            throw matrix_error();
        }
        return data[i - 1][j];
    }

    // TODO: 重载乘号，返回矩阵乘法 lhs * rhs 的结果。如果 lhs 的列数与 rhs 的行数不相等，抛出 matrix_error 错误。
    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m) {
            throw matrix_error();
        }
        matrix res(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; ++i) {
            for (int j = 0; j < rhs.n; ++j) {
                fraction sum(0);
                for (int k = 0; k < lhs.n; ++k) {
                    sum = sum + lhs.data[i][k] * rhs.data[k][j];
                }
                res.data[i][j] = sum;
            }
        }
        return res;
    }

    // TODO: 返回矩阵的转置。若矩阵为空，抛出 matrix_error 错误。
    matrix transposition() {
        if (m == 0 || n == 0) {
            throw matrix_error();
        }
        matrix res(n, m);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                res.data[j][i] = data[i][j];
            }
        }
        return res;
    }

    // TODO: 返回矩阵的行列式。建议用高斯消元实现。若矩阵不是方阵或为空，抛出 matrix_error 错误。
    fraction determination() {
        if (m != n || m == 0) {
            throw matrix_error();
        }
        matrix tmp(*this);
        fraction det(1);
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i; j < n; ++j) {
                if (!(tmp.data[j][i] == fraction(0))) {
                    pivot = j;
                    break;
                }
            }
            if (tmp.data[pivot][i] == fraction(0)) {
                return fraction(0);
            }
            if (pivot != i) {
                std::swap(tmp.data[i], tmp.data[pivot]);
                det = det * fraction(-1);
            }
            det = det * tmp.data[i][i];
            fraction inv = fraction(1) / tmp.data[i][i];
            for (int j = i + 1; j < n; ++j) {
                fraction factor = tmp.data[j][i] * inv;
                for (int k = i; k < n; ++k) {
                    tmp.data[j][k] = tmp.data[j][k] - factor * tmp.data[i][k];
                }
            }
        }
        return det;
    }
};

class resistive_network {
private:

    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 矩阵A 和 矩阵C
    matrix adjacency, conduction;

    // 预计算的 A^T C A
    matrix AT_C_A;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************
    matrix get_submatrix(const matrix& mat, const std::vector<int>& rows_to_remove, const std::vector<int>& cols_to_remove) {
        int new_m = mat.m - rows_to_remove.size();
        int new_n = mat.n - cols_to_remove.size();
        matrix res(new_m, new_n);
        int r = 0;
        for (int i = 0; i < mat.m; ++i) {
            bool remove_r = false;
            for (int rr : rows_to_remove) {
                if (i == rr) { remove_r = true; break; }
            }
            if (remove_r) continue;
            int c = 0;
            for (int j = 0; j < mat.n; ++j) {
                bool remove_c = false;
                for (int cc : cols_to_remove) {
                    if (j == cc) { remove_c = true; break; }
                }
                if (remove_c) continue;
                res.data[r][c] = mat.data[i][j];
                c++;
            }
            r++;
        }
        return res;
    }

public:

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

    // TODO: 设置电阻网络，构建矩阵A和C。节点数量为interface_size_，接线数量为connection_size_。
    //       对于 1<=i<=connection_size_，从节点from[i-1]到节点to[i-1]有接线，对应电阻为resistance[i-1]。
    //       保证接线使得电阻网络联通，from[i-1] < to[i-1]，resitance[i-1] > 0，均合法。
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;
        adjacency = matrix(connection_size, interface_size);
        conduction = matrix(connection_size, connection_size);
        for (int i = 0; i < connection_size; ++i) {
            adjacency(i + 1, from[i] - 1) = fraction(1);
            adjacency(i + 1, to[i] - 1) = fraction(-1);
            conduction(i + 1, i) = fraction(1) / resistance[i];
        }
        AT_C_A = adjacency.transposition() * conduction * adjacency;
    }

    ~resistive_network() = default;

    // TODO: 返回节点 interface_id1 和 interface_id2 (1-based)之间的等效电阻。
    //       保证 interface_id1 <= interface_id2 均合法。
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) return fraction(0);
        
        std::vector<int> remove_i = {interface_id1 - 1};
        matrix Mi = get_submatrix(AT_C_A, remove_i, remove_i);
        
        std::vector<int> remove_ij = {interface_id1 - 1, interface_id2 - 1};
        matrix Mij = get_submatrix(AT_C_A, remove_ij, remove_ij);
        
        fraction det_Mij = (Mij.m == 0) ? fraction(1) : Mij.determination();
        return det_Mij / Mi.determination();
    }

    // TODO: 在给定节点电流I的前提下，返回节点id(1-based)的电压。认为节点interface_size(1-based)的电压为0。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电流为 current[i-1]。
    //       保证 current 使得电阻网络有解，id < interface_size 合法。
    fraction get_voltage(int id, fraction current[]) {
        std::vector<int> remove_n = {interface_size - 1};
        matrix Mn = get_submatrix(AT_C_A, remove_n, remove_n);
        
        matrix AT_C_A_replaced = AT_C_A;
        for (int i = 0; i < interface_size; ++i) {
            AT_C_A_replaced.data[i][id - 1] = current[i];
        }
        matrix Mni = get_submatrix(AT_C_A_replaced, remove_n, remove_n);
        
        return Mni.determination() / Mn.determination();
    }


    // TODO: 在给定节点电压U的前提下，返回电阻网络的功率。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电压为 voltage[i-1]。
    //       保证 voltage 合法。
    fraction get_power(fraction voltage[]) {
        matrix U(interface_size, 1);
        for (int i = 0; i < interface_size; ++i) {
            U.data[i][0] = voltage[i];
        }
        matrix uw = adjacency * U;
        fraction P(0);
        for (int i = 0; i < connection_size; ++i) {
            fraction uwi = uw.data[i][0];
            P = P + uwi * uwi * conduction.data[i][i];
        }
        return P;
    }
};


#endif //SRC_HPP