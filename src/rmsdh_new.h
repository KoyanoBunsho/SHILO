#ifndef rmsdh_new
#define rmsdh_new
#include "rmsd_io.h"
#include "rmsd_struct.h"
#include <chrono>
#include <float.h>
#include <functional>
#include <math.h>
#include <memory>
#include <numeric>
#include <omp.h>
#include <unordered_map>
constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double eps = std::numeric_limits<double>::epsilon();
// https://sotanishy.github.io/cp-library-cpp/dp/smawk.cpp
template <typename T>
std::vector<int> smawk(int n, int m, std::function<T(int, int)> lookUpij) {
  auto smawkSub = [&](auto &smawkSub, const std::vector<int> &row,
                      const std::vector<int> &col) -> std::vector<int> {
    int row_size = row.size();
    if (row_size == 0) {
      return {};
    }
    std::vector<int> col2;
    for (int c : col) {
      while (!col2.empty() && lookUpij(row[col2.size() - 1], col2.back()) >
                                  lookUpij(row[col2.size() - 1], c)) {
        col2.pop_back();
      }
      if ((int)col2.size() < row_size) {
        col2.push_back(c);
      }
    }
    std::vector<int> row2;
    for (int i = 1; i < row_size; i += 2) {
      row2.push_back(row[i]);
    }
    auto sub = smawkSub(smawkSub, row2, col2);
    std::vector<int> ans(row_size);
    int sub_size = (int)sub.size();
    for (int i = 0; i < sub_size; ++i) {
      ans[2 * i + 1] = sub[i];
    }
    int j = 0;
    for (int i = 0; i < row_size; i += 2) {
      ans[i] = col2[j];
      int end = (i == row_size - 1 ? col2.back() : ans[i + 1]);
      while (col2[j] < end) {
        ++j;
        if (lookUpij(row[i], ans[i]) > lookUpij(row[i], col2[j])) {
          ans[i] = col2[j];
        }
      }
    }
    return ans;
  };

  std::vector<int> row(n), col(m);
  std::iota(row.begin(), row.end(), 0);
  std::iota(col.begin(), col.end(), 0);
  return smawkSub(smawkSub, row, col);
}

class ProteinRMSDhinge {
public:
  ProteinRMSDhinge(const Eigen::Matrix3Xd &A, const Eigen::Matrix3Xd &B,
                   const int hinge_num) {
    assert(A.cols() == B.cols());
    int n = A.cols();
    this->A = A;
    this->B = B;
    this->n = n;
    this->hinge_num = hinge_num;
    this->aTa_acc = std::vector<double>(n + 1);
    this->bTb_acc = std::vector<double>(n + 1);
    this->baT_acc = std::vector<Eigen::Matrix3Xd>(n + 1);
    this->a_acc = std::vector<Eigen::Vector3d>(n + 1);
    this->b_acc = std::vector<Eigen::Vector3d>(n + 1);
    this->a_acc[0] = Eigen::Vector3d::Zero(3);
    this->b_acc[0] = Eigen::Vector3d::Zero(3);
    this->aTa_acc[0] = 0;
    this->bTb_acc[0] = 0;
    this->baT_acc[0] = Eigen::MatrixXd::Zero(3, 3);
    this->A_T = A.transpose();
    this->B_T = B.transpose();
    for (int i = 0; i < n; i++) {
      this->a_acc[i + 1] = this->a_acc[i] + this->A.col(i);
      this->b_acc[i + 1] = this->b_acc[i] + this->B.col(i);
      this->aTa_acc[i + 1] =
          this->aTa_acc[i] + this->A_T.row(i) * this->A.col(i);
      this->bTb_acc[i + 1] =
          this->bTb_acc[i] + this->B_T.row(i) * this->B.col(i);
      this->baT_acc[i + 1] =
          this->baT_acc[i] + this->B.col(i) * this->A_T.row(i);
    }
  }
  bool left_is_leq_right(double left, double right) {
    return left - right <= DBL_EPSILON * fmax(1, fmax(fabs(left), fabs(right)));
  }
  inline Eigen::Vector3d
  calcCentroid(const std::vector<Eigen::Vector3d> &val_acc, int i, int j) {
    /*
    Calculate the centroid of a 3d structure in constant time
    val_acc: Accumulative sum of 3d vector
    val_acc[j] - val_acc[i - 1] represents the sum of 3d vecor from the i to j
    Centroid is the average coordinate of all coordinates betwee i and j
    */
    Eigen::Vector3d val_bar = (val_acc[j] - val_acc[i - 1]) / (j - i + 1);
    return val_bar;
  }

  inline double calcSubSD(const double &pTp_sub, const double &qTq_sub,
                          const Eigen::Matrix3Xd &H_sub) {
    /*
    Calculate SD of a coordination pair between the index i and j
    sd = p^Tp + q^Tq - 2 trace(RH)
    */
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H_sub, Eigen::ComputeThinU |
                                                     Eigen::ComputeThinV);
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    Eigen::Matrix3d V = svd.matrixV();
    Eigen::Matrix3d U_T = svd.matrixU().transpose();
    double d = (V * U_T).determinant();
    (d > 0.0) ? d = 1.0 : d = -1.0;
    I(2, 2) = d;
    Eigen::MatrixXd RH = V * I * U_T * H_sub; // R = V * I * U_T
    double sub_sd = pTp_sub + qTq_sub - 2 * RH.trace();
    (sub_sd > 0.0) ? sub_sd = sub_sd : sub_sd = 0.0;
    return sub_sd;
  }

  inline double calcSDij(int i, int j) {
    double aTa_sub = aTa_acc[j] - aTa_acc[i - 1];
    double bTb_sub = bTb_acc[j] - bTb_acc[i - 1];
    Eigen::Matrix3Xd baT_sub = baT_acc[j] - baT_acc[i - 1];
    Eigen::Vector3d a_sub = a_acc[j] - a_acc[i - 1];
    Eigen::Vector3d b_sub = b_acc[j] - b_acc[i - 1];
    Eigen::Vector3d a_bar = calcCentroid(a_acc, i, j);
    Eigen::Vector3d b_bar = calcCentroid(b_acc, i, j);
    Eigen::Matrix3Xd H_sub = baT_sub - b_sub * a_bar.transpose() -
                             b_bar * a_sub.transpose() +
                             (j - i + 1) * b_bar * a_bar.transpose();
    double pTp_sub = aTa_sub - a_sub.transpose() * a_bar -
                     a_bar.transpose() * a_sub +
                     (j - i + 1) * a_bar.transpose() * a_bar;
    double qTq_sub = bTb_sub - b_sub.transpose() * b_bar -
                     b_bar.transpose() * b_sub +
                     (j - i + 1) * b_bar.transpose() * b_bar;
    double sub_sd_ij = calcSubSD(pTp_sub, qTq_sub, H_sub);
    return sub_sd_ij;
  }
  void initializeDPandBackTrackingTable(
      std::vector<std::vector<double>> &dp,
      std::vector<std::vector<int>> &backtracking_table) {
    for (int j = 0; j < n + 1; j++) {
      dp[j][1] = INF;
    }
    backtracking_table[1][1] = 0;
    for (int j = 1; j < n; j++) {
      double sub_sd_1j = calcSDij(1, j);
      dp[j][0] = sub_sd_1j;
    }
  }
  std::vector<int>
  backTrackDPKFor(std::vector<std::vector<int>> backtracking_table,
                  std::vector<int> &hinge_index_vec, int j_index) {
    for (int k = hinge_num; k > 0; k--) {
      j_index = backtracking_table[j_index - 1][k];
      hinge_index_vec.push_back(j_index);
    }
    return hinge_index_vec;
  }

  RMSDhHingeCnt CalcRMSDhK() {
    /*
    k (the number of hinges) is expected to be given in this function
    A: Protein coordinates matrix (3 x n)
    B: Protein coordinates matrix (3 x n)
    n: the number of coordinate pairs
    hinge_num: hinge number to find

    DP definition is as follows.
    OPT(j, k) = min_{1 <= i <= j} {sd_{ij} + OPT(i - 1, k-1)}
    Final RMSDh = sqrt((OPT(j, k)) / n)
    k: the number of hinges

    Precomputed values are as follows
    a, b, a^T * a, b^T * b, b * a^T

    H    = qp^T
         = (b - b_bar)(a - a_bar)^T
         = b * a^T - b * a_bar^T - b_bar * a^T + b_bar * a_bar^T
    p^Tp = (a - a_bar)^T(a - a_bar)
         =  a^T * a - a^T * a_bar - a_bar^T * a + a_bar^T * a_bar
    q^Tq = (b - b_bar)^T(b - b_bar)
         =  b^T * b - b^T * b_bar - b_bar^T * b + b_bar^T * b_bar
    */
    std::vector<std::vector<double>> dp(n + 1,
                                        std::vector<double>(hinge_num + 1));
    std::vector<int> hinge_index_vec;
    std::vector<std::vector<int>> backtracking_table(
        n + 1, std::vector<int>(hinge_num + 1));
    initializeDPandBackTrackingTable(dp, backtracking_table);
    for (int k = 1; k < hinge_num + 1; k++) {
      for (int j = 1; j < n + 1; j++) {
        double ans = INF;
        int min_index = 1;
        for (int i = k; i <= j; i++) {
          double sub_sd = calcSDij(i, j);
          // if (left_is_leq_right(sub_sd + dp[i - 1][k - 1], ans)) {
          if (sub_sd + dp[i - 1][k - 1] < ans) {
            ans = sub_sd + dp[i - 1][k - 1];
            min_index = i;
          }
        }
        dp[j][k] = ans;
        backtracking_table[j][k] = min_index;
      }
    }
    double rmsdh_result;
    hinge_index_vec =
        backTrackDPKFor(backtracking_table, hinge_index_vec, n + 1);
    RMSDhHingeCnt rmsdh_hinge_cnt_result;
    if (dp[n][hinge_num] < 0) {
      rmsdh_result = 0.0;
    } else {
      rmsdh_result = std::sqrt(dp[n][hinge_num] / n);
    }
    rmsdh_hinge_cnt_result.rmsdh_result = rmsdh_result;
    rmsdh_hinge_cnt_result.hinge_cnt = hinge_num;
    rmsdh_hinge_cnt_result.hinge_index_vec = hinge_index_vec;
    return rmsdh_hinge_cnt_result;
  }
  void updateHingePosition(std::vector<int> &hinge_index_vec,
                           const int hinge_num, int hinge_start) {
    for (int i = hinge_start; i < hinge_num + 1; i += 2) {
      int hinge_left, hinge_right;
      if (i == 1) {
        hinge_left = 1;
      } else {
        hinge_left = hinge_index_vec[i - 2];
      }
      if (i == hinge_num) {
        hinge_right = n + 1;
      } else {
        hinge_right = hinge_index_vec[i];
      }
      int corrected_hinge_index = hinge_index_vec[i - 1];
      double ans = calcSDij(hinge_left, corrected_hinge_index - 1) +
                   calcSDij(corrected_hinge_index, hinge_right - 1);
      for (int hinge_index = hinge_left; hinge_index < hinge_right - 1;
           hinge_index++) {
        double sub_sd_left = calcSDij(hinge_left, hinge_index);
        double sub_sd_right = calcSDij(hinge_index + 1, hinge_right - 1);
        if (ans > sub_sd_left + sub_sd_right) {
          ans = sub_sd_left + sub_sd_right;
          corrected_hinge_index = hinge_index + 1;
        }
      }
      hinge_index_vec[i - 1] = corrected_hinge_index;
    }
  }

  RMSDhHingeCnt calcRMSDhKAfterHingeUpdate(std::vector<int> &hinge_index_vec,
                                           const int hinge_num) {
    double rmsdh_result = 0.0;
    std::sort(hinge_index_vec.begin(), hinge_index_vec.end());
    hinge_index_vec.push_back(n + 1);
    hinge_index_vec.insert(hinge_index_vec.begin(), 1);
    for (int i = 0; i < (int)hinge_index_vec.size() - 1; i++) {
      double sub_sd = calcSDij(hinge_index_vec[i], hinge_index_vec[i + 1] - 1);
      rmsdh_result += sub_sd;
    }
    rmsdh_result = std::sqrt(rmsdh_result / static_cast<double>(n));
    RMSDhHingeCnt rmsdh_hinge_cnt_result;
    hinge_index_vec.pop_back();
    hinge_index_vec.erase(hinge_index_vec.begin());
    rmsdh_hinge_cnt_result.rmsdh_result = rmsdh_result;
    rmsdh_hinge_cnt_result.hinge_cnt = hinge_num;
    rmsdh_hinge_cnt_result.hinge_index_vec = hinge_index_vec;
    return rmsdh_hinge_cnt_result;
  }
  bool updateHingePositionBackward(std::vector<int> &hinge_index_vec,
                                   const int hinge_num) {
    bool is_updated = false;
    for (int i = 0; i < hinge_num; i++) {
      int hinge_left;
      if (i == hinge_num - 1) {
        hinge_left = 1;
      } else {
        hinge_left = hinge_index_vec[hinge_num - i - 2];
      }
      int hinge_right;
      if (i == 0) {
        hinge_right = n + 1;
      } else {
        hinge_right = hinge_index_vec[hinge_num - i];
      }
      int corrected_hinge_index = hinge_index_vec[hinge_num - i - 1];
      double ans = calcSDij(hinge_left, corrected_hinge_index - 1) +
                   calcSDij(corrected_hinge_index, hinge_right - 1);
      for (int hinge_index = hinge_left; hinge_index < hinge_right - 1;
           hinge_index++) {
        double sub_sd_left = calcSDij(hinge_left, hinge_index);
        double sub_sd_right = calcSDij(hinge_index + 1, hinge_right - 1);
        if (ans > sub_sd_left + sub_sd_right) {
          ans = sub_sd_left + sub_sd_right;
          corrected_hinge_index = hinge_index + 1;
          is_updated = true;
        }
      }
      hinge_index_vec[hinge_num - i - 1] = corrected_hinge_index;
    }
    return is_updated;
  }
  void updateHingePositionForward(std::vector<int> &hinge_index_vec,
                                  const int hinge_num) {
    for (int i = 1; i < hinge_num; i++) {
      int hinge_left = hinge_index_vec[i - 1];
      int hinge_right;
      if (i == hinge_num - 1) {
        hinge_right = n + 1;
      } else {
        hinge_right = hinge_index_vec[i + 1];
      }
      int corrected_hinge_index = hinge_index_vec[i];
      double ans = calcSDij(hinge_left, corrected_hinge_index - 1) +
                   calcSDij(corrected_hinge_index, hinge_right - 1);
      for (int hinge_index = hinge_left; hinge_index < hinge_right - 1;
           hinge_index++) {
        double sub_sd_left = calcSDij(hinge_left, hinge_index);
        double sub_sd_right = calcSDij(hinge_index + 1, hinge_right - 1);
        if (ans > sub_sd_left + sub_sd_right) {
          ans = sub_sd_left + sub_sd_right;
          corrected_hinge_index = hinge_index + 1;
        }
      }
      hinge_index_vec[i] = corrected_hinge_index;
    }
  }

  RMSDhHingeCnt RMSDhkPostProcessing(std::vector<int> &hinge_index_vec,
                                     const int hinge_num) {
    std::sort(hinge_index_vec.begin(), hinge_index_vec.end());
    bool is_updated = updateHingePositionBackward(hinge_index_vec, hinge_num);
    if (is_updated) {
      updateHingePositionForward(hinge_index_vec, hinge_num);
    }
    return calcRMSDhKAfterHingeUpdate(hinge_index_vec, hinge_num);
  }
  AblationResult RMSDhkPostProcessingLoop(std::vector<int> &hinge_index_vec,
                                          const int hinge_num) {
    std::sort(hinge_index_vec.begin(), hinge_index_vec.end());
    bool is_updated = true;
    int updated_cnt = 0;
    int max_cnt = 10;
    while (is_updated && updated_cnt <= max_cnt) {
      updated_cnt++;
      is_updated = updateHingePositionBackward(hinge_index_vec, hinge_num);
      if (is_updated) {
        updateHingePositionForward(hinge_index_vec, hinge_num);
      }
    }
    AblationResult res;
    res.iter_num = updated_cnt;
    RMSDhHingeCnt rmsdhk =
        calcRMSDhKAfterHingeUpdate(hinge_index_vec, hinge_num);
    res.rmsdh_result = rmsdhk.rmsdh_result;
    res.hinge_cnt = rmsdhk.hinge_cnt;
    res.hinge_index_vec = rmsdhk.hinge_index_vec;
    return res;
  }

  bool MonotonicityChecker(double a, double b, double c, double d) {
    if (a == b && b == c && c == d) {
      return true;
    }
    if (c <= d) {
      if (a <= b) {
        return true;
      } else {
        return false;
      }
    }
    if (b < a) {
      if (d < c) {
        return true;
      } else {
        return false;
      }
    } else {
      return true;
    }
  }

  double calcMongeRate() {
    int monge_cnt = 0;
    int total_cnt = 0;
    int monotonicity_cnt = 0;
    for (int j = 1; j < n; j++) {
      for (int i = 1; i <= j; i++) {
        total_cnt++;
        if (calcSDij(i, j) + calcSDij(i + 1, j + 1) <=
            calcSDij(i + 1, j) + calcSDij(i, j + 1)) {
          monge_cnt++;
        }
        if (MonotonicityChecker(calcSDij(i, j), calcSDij(i, j + 1),
                                calcSDij(i + 1, j), calcSDij(i + 1, j + 1))) {
          monotonicity_cnt++;
        }
      }
    }
    std::cout << static_cast<double>(monotonicity_cnt) /
                     static_cast<double>(total_cnt)
              << std::endl;
    assert(total_cnt == n * (n - 1) / 2);
    return static_cast<double>(monge_cnt) / static_cast<double>(total_cnt);
  }
  double calcMongeRate2() {
    std::vector<std::vector<double>> tmr_dp(n + 1, std::vector<double>(n + 1));
    for (int i = 1; i <= n; i++) {
      for (int j = i; j <= n; j++) {
        tmr_dp[j][i] = calcSDij(i, j);
      }
    }
    long long monotonicity_cnt = 0;
    long long total_cnt = 0;
#pragma omp parallel for reduction(+:monotonicity_cnt, total_cnt) num_threads(1000)
    for (int j = 2; j <= n; j++) {
      for (int i = 1; i < j; i++) {
        for (int l = 2; l < i; l++) {
          for (int k = 1; k < l; k++) {
            total_cnt++;
            assert(k < l && l < i && i < j);
            if (tmr_dp[i][k] + tmr_dp[j][l] <= tmr_dp[i][l] + tmr_dp[j][k]) {
              monotonicity_cnt++;
            }
          }
        }
      }
    }
    // long long lln = static_cast<long long>(n);
    // assert(total_cnt ==
    //     lln * (lln - 1) * (lln - 2) * (lln - 3) / (4 * 3 * 2));
    std::cout << "total_cnt: " << total_cnt
              << ", monotonicity_cnt: " << monotonicity_cnt << std::endl;
    double tmr =
        static_cast<double>(monotonicity_cnt) / static_cast<double>(total_cnt);
    return tmr;
  }
  double calcMongeRateFromMatrix(std::vector<std::vector<double>> &tmr_dp) {
    int monge_cnt = 0;
    int total_cnt = 0;
    for (int j = 1; j < n; j++) {
      for (int i = 1; i <= j; i++) {
        total_cnt++;
        if (tmr_dp[j][i] + tmr_dp[j + 1][i + 1] <=
            tmr_dp[j + 1][i] + tmr_dp[j][i + 1]) {
          monge_cnt++;
        }
      }
    }
    assert(total_cnt == n * (n - 1) / 2);
    return static_cast<double>(monge_cnt) / static_cast<double>(total_cnt);
  }
  RMSDhHingeCnt CalcFastRMSDhK(bool is_postprocessing = false) {
    /*
    Add SMAWK algorithm to compute RMSDhk in O(kn)
    k (the number of hinges) is expected to be given in this function
    A: Protein coordinates matrix (3 x n)
    B: Protein coordinates matrix (3 x n)
    n: the number of coordinate pairs
    hinge_num: hinge number to find

    DP definition is as follows.
    OPT(j, k) = min_{1 <= i <= j} {sd_{ij} + OPT(i - 1, k-1)}
    Final RMSDh = sqrt((OPT(j, k)) / n)
    k: the number of hinges

    Precomputed values are as follows
    a, b, a^T * a, b^T * b, b * a^T

    H    = qp^T
         = (b - b_bar)(a - a_bar)^T
         = b * a^T - b * a_bar^T - b_bar * a^T + b_bar * a_bar^T
    p^Tp = (a - a_bar)^T(a - a_bar)
         =  a^T * a - a^T * a_bar - a_bar^T * a + a_bar^T * a_bar
    q^Tq = (b - b_bar)^T(b - b_bar)
         =  b^T * b - b^T * b_bar - b_bar^T * b + b_bar^T * b_bar
    */
    std::vector<std::vector<double>> dp(n + 1,
                                        std::vector<double>(hinge_num + 1));
    std::vector<int> hinge_index_vec;
    std::vector<std::vector<int>> backtracking_table(
        n + 1, std::vector<int>(hinge_num + 1));
    initializeDPandBackTrackingTable(dp, backtracking_table);
    for (int k = 1; k < hinge_num + 1; k++) {
      auto lookUpij = [&](int j, int i) -> double {
        // change 0-index to 1-index
        // i and j are 0-index.
        if (i > j) {
          return INF;
        } else {
          double sub_sd = calcSDij(i + 1, j + 1);
          return sub_sd + dp[i][k - 1];
        }
      };
      std::vector<int> argmin = smawk<double>(n, n, lookUpij);
      int argmin_size = (int)argmin.size();
      for (int j = 0; j < argmin_size; j++) {
        // argmin[j] and j are 0-index, and the others are 1-index
        backtracking_table[j + 1][k] = argmin[j] + 1;
        double sub_sd = calcSDij(argmin[j] + 1, j + 1);
        dp[j + 1][k] = dp[argmin[j]][k - 1] + sub_sd;
      }
    }
    double rmsdh_result;
    hinge_index_vec =
        backTrackDPKFor(backtracking_table, hinge_index_vec, n + 1);
    int hinge_index_vec_size = (int)hinge_index_vec.size();
    for (int i = 0; i < hinge_index_vec_size; i++) {
      hinge_index_vec[i];
    }
    RMSDhHingeCnt rmsdh_hinge_cnt_result;
    if (is_postprocessing) {
      RMSDhHingeCnt rmsdh_hinge_cnt_result =
          RMSDhkPostProcessing(hinge_index_vec, hinge_num);
      return rmsdh_hinge_cnt_result;
    } else {
      rmsdh_hinge_cnt_result.hinge_index_vec = hinge_index_vec;
    }
    if (dp[n][hinge_num] < 0) {
      rmsdh_result = 0.0;
    } else {
      rmsdh_result = std::sqrt(dp[n][hinge_num] / (double)n);
    }
    rmsdh_hinge_cnt_result.rmsdh_result = rmsdh_result;
    rmsdh_hinge_cnt_result.hinge_cnt = hinge_num;
    return rmsdh_hinge_cnt_result;
  }
  AblationResult CalcFastRMSDhKLoop() {
    /*
    Add SMAWK algorithm to compute RMSDhk in O(kn)
    k (the number of hinges) is expected to be given in this function
    A: Protein coordinates matrix (3 x n)
    B: Protein coordinates matrix (3 x n)
    n: the number of coordinate pairs
    hinge_num: hinge number to find

    DP definition is as follows.
    OPT(j, k) = min_{1 <= i <= j} {sd_{ij} + OPT(i - 1, k-1)}
    Final RMSDh = sqrt((OPT(j, k)) / n)
    k: the number of hinges

    Precomputed values are as follows
    a, b, a^T * a, b^T * b, b * a^T

    H    = qp^T
         = (b - b_bar)(a - a_bar)^T
         = b * a^T - b * a_bar^T - b_bar * a^T + b_bar * a_bar^T
    p^Tp = (a - a_bar)^T(a - a_bar)
         =  a^T * a - a^T * a_bar - a_bar^T * a + a_bar^T * a_bar
    q^Tq = (b - b_bar)^T(b - b_bar)
         =  b^T * b - b^T * b_bar - b_bar^T * b + b_bar^T * b_bar
    */
    std::vector<std::vector<double>> dp(n + 1,
                                        std::vector<double>(hinge_num + 1));
    std::vector<int> hinge_index_vec;
    std::vector<std::vector<int>> backtracking_table(
        n + 1, std::vector<int>(hinge_num + 1));
    initializeDPandBackTrackingTable(dp, backtracking_table);
    for (int k = 1; k < hinge_num + 1; k++) {
      auto lookUpij = [&](int j, int i) -> double {
        // change 0-index to 1-index
        // i and j are 0-index.
        if (i > j) {
          return INF;
        } else {
          double sub_sd = calcSDij(i + 1, j + 1);
          return sub_sd + dp[i][k - 1];
        }
      };
      std::vector<int> argmin = smawk<double>(n, n, lookUpij);
      int argmin_size = (int)argmin.size();
      for (int j = 0; j < argmin_size; j++) {
        // argmin[j] and j are 0-index, and the others are 1-index
        backtracking_table[j + 1][k] = argmin[j] + 1;
        double sub_sd = calcSDij(argmin[j] + 1, j + 1);
        dp[j + 1][k] = dp[argmin[j]][k - 1] + sub_sd;
      }
    }
    hinge_index_vec =
        backTrackDPKFor(backtracking_table, hinge_index_vec, n + 1);
    int hinge_index_vec_size = (int)hinge_index_vec.size();
    AblationResult rmsdh_hinge_cnt_result;
    rmsdh_hinge_cnt_result =
        RMSDhkPostProcessingLoop(hinge_index_vec, hinge_num);
    return rmsdh_hinge_cnt_result;
  }

private:
  Eigen::Matrix3Xd A, B;
  int n;
  std::vector<double> aTa_acc, bTb_acc;
  std::vector<Eigen::Matrix3Xd> baT_acc;
  std::vector<Eigen::Vector3d> a_acc, b_acc;
  Eigen::MatrixXd A_T, B_T;
  int hinge_num;
};

double CalcSD(Eigen::Matrix3Xd P, Eigen::Matrix3Xd Q, std::vector<double> W) {
  if (P.cols() != Q.cols()) {
    throw "CalcRMSD(): input data mis-match";
  }
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 3);
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  int p_cols_num = P.cols();
  for (int i = 0; i < p_cols_num; i++) {
    H += W[i] * Q.col(i) * P.col(i).transpose();
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU |
                                               Eigen::ComputeThinV);
  Eigen::Matrix3d V = svd.matrixV();
  Eigen::Matrix3d U_T = svd.matrixU().transpose();
  double d = (V * U_T).determinant();
  (d > 0.0) ? d = 1.0 : d = -1.0;
  I(2, 2) = d;
  Eigen::Matrix3d R = V * I * U_T;
  double sd = 0;
  for (int i = 0; i < p_cols_num; i++) {
    sd += W[i] * P.col(i).transpose() * P.col(i);
    sd += W[i] * Q.col(i).transpose() * Q.col(i);
  }
  Eigen::MatrixXd RH = R * H;
  sd -= 2 * RH.trace();
  if (sd < 0) {
    sd = 0;
  }
  return sd;
}

double CalcRMSD(Eigen::Matrix3Xd P, Eigen::Matrix3Xd Q, std::vector<double> W) {
  double sd = CalcSD(P, Q, W);
  double rmsd_result;
  int p_cols_num = P.cols();
  rmsd_result = std::sqrt(sd / p_cols_num);
  return rmsd_result;
}

ConformationPair MoveToOrigin(Eigen::Matrix3Xd P, Eigen::Matrix3Xd Q,
                              std::vector<double> W) {
  Eigen::Vector3d p{0, 0, 0}, q{0, 0, 0};
  ConformationPair PQ_pair;
  double w_sum = 0.0;
  int p_cols_num = P.cols();
  for (int i = 0; i < p_cols_num; i++) {
    w_sum += W[i];
  }
  for (int i = 0; i < p_cols_num; i++) {
    p += (W[i] * P.col(i)) / w_sum;
  }
  for (int i = 0; i < p_cols_num; i++) {
    q += (W[i] * Q.col(i)) / w_sum;
  }
  Eigen::Matrix3Xd X = P.colwise() - p;
  Eigen::Matrix3Xd Y = Q.colwise() - q;
  PQ_pair.P = X;
  PQ_pair.Q = Y;
  return PQ_pair;
}

#endif
