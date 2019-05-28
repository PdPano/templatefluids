#include "../../grid/cartesian_grid.hpp"
#include "../filters/minimal_filter_factory.hpp"
#include "../flag_handler.hpp"
#include "luisa_detector_345.hpp"
#include <iostream>

LuisaDetector345::LuisaDetector345(double sensitivity_in, int nPointsI_in,
                                   int nPointsJ_in)
    : sensitivity(sensitivity_in), nPointsI(nPointsI_in), nPointsJ(nPointsJ_in),
      rho_to_derive(std::vector<double>(std::max(nPointsI, nPointsJ))),
      rho_to_derive_two(std::vector<double>(std::max(nPointsI, nPointsJ))),
      rho_to_load(std::vector<double>(std::max(nPointsI, nPointsJ))),
      numerator(std::vector<double>(std::max(nPointsI, nPointsJ))),
      denominator(std::vector<double>(std::max(nPointsI, nPointsJ))),
      der_flags(std::vector<int>(std::max(nPointsI, nPointsJ))),
      minimal_filter(create_minimal_filter(3)) {}

void LuisaDetector345::detect_shocks(const CartesianGrid &grid) {
  shocked_points.clear();
  for (int line = 0; line < nPointsI; line++) {
    load_data(line, 'x', grid);
    minimal_filter->filter_line(rho_to_load, &rho_to_derive, der_flags,
                                nPointsJ, 1);
    minimal_filter->filter_line(rho_to_load, &rho_to_derive_two, der_flags,
                                nPointsJ, 1);
    find_shocks(line, 'x', &shocked_points, grid);
  }
  for (int col = 0; col < nPointsJ; col++) {
    load_data(col, 'y', grid);
    minimal_filter->filter_line(rho_to_load, &rho_to_derive, der_flags,
                                nPointsI, 1);
    minimal_filter->filter_line(rho_to_load, &rho_to_derive_two, der_flags,
                                nPointsI, 1);
    find_shocks(col, 'y', &shocked_points, grid);
  }
}

void LuisaDetector345::load_data(int line_or_col, char dir,
                                 const CartesianGrid &grid) {
  int len;
  int ind;
  int shift;

  if (dir == 'x') {
    len = nPointsJ;
    ind = grid.IND(line_or_col, 0);
    shift = grid.shiftX();
  } else {
    len = nPointsI;
    ind = grid.IND(0, line_or_col);
    shift = grid.shiftY();
  }

  double max_rho = -1e20;
  double min_rho = 1e20;
  for (int i = 0; i < len; i++) {
    double rho = grid.rho(ind);
    if (min_rho > rho) {
      min_rho = rho;
    }
    if (max_rho < rho) {
      max_rho = rho;
    }
    rho_to_load[i] = rho;
    if (dir == 'x') {
      der_flags[i] = flag_functions::derx(grid.flag(ind));
    } else {
      der_flags[i] = flag_functions::dery(grid.flag(ind));
    }
    ind += shift;
  }
  if (max_rho - min_rho < 1e-2) {
    max_rho = min_rho + 1e-2;
  }
  for (int i = 0; i < len; i++) {
    rho_to_load[i] = (rho_to_load[i] - min_rho) / (max_rho - min_rho) +
                     0.1 * sin(i * M_PI / len);
  }
}

void LuisaDetector345::find_shocks(int line_or_col, char dir,
                                   std::set<int> *shocked_points,
                                   const CartesianGrid &grid) {
  int len;
  if (dir == 'x') {
    len = nPointsJ;
  } else {
    len = nPointsI;
  }

  auto convert_ind = [&](int i) {
    if (dir == 'x') {
      return grid.IND(line_or_col, i);
    }
    return grid.IND(i, line_or_col);
  };

  compute_numerator(len);
  compute_denominator(len);

  for (int i = 0; i < len; i++) {
    auto ratio = (numerator[i] + 1e-80) / (denominator[i] + 1e-80);
    auto continuous_up_to = -log(ratio) / log(2.0);
    auto l_flag = flag_functions::left(der_flags[i]);
    auto r_flag = flag_functions::right(der_flags[i]);
    if (continuous_up_to < sensitivity or std::min(l_flag, r_flag) < 5) {
      shocked_points->insert(convert_ind(i));
    }
  }
  std::set<int> expand;
  for (int i = 2; i < len - 2; i++) {
    auto shock_exists = [&](int ind) {
      return shocked_points->find(ind) != shocked_points->end();
    };
    if (shock_exists(convert_ind(i - 2)) or shock_exists(convert_ind(i - 1)) or
        shock_exists(convert_ind(i + 1)) or shock_exists(convert_ind(i + 2))) {
      expand.insert(convert_ind(i));
    }
  }
  for (auto &i : expand) {
    shocked_points->insert(i);
  }
}

void LuisaDetector345::derive_vector(std::vector<double> *in, int len,
                                     int step) {
  auto &aux = rho_to_load;
  for (int i = 0; i < len; i++) {
    aux[i] = generic_first_der(i, *in, step);
  }
  in->swap(aux);
}

void LuisaDetector345::compute_numerator(int len) {
  derive_vector(&rho_to_derive, len, 1); // drho
  derive_vector(&rho_to_derive, len, 1); // d2rho
  derive_vector(&rho_to_derive, len, 1); // d3rho
  for (int i = 0; i < len; i++) {
    numerator[i] = fabs(rho_to_derive[i]);
  }
  derive_vector(&rho_to_derive, len, 1); // d4rho
  for (int i = 0; i < len; i++) {
    numerator[i] += fabs(rho_to_derive[i]);
  }
  derive_vector(&rho_to_derive, len, 1); // d5rho
  for (int i = 0; i < len; i++) {
    numerator[i] += fabs(rho_to_derive[i]);
  }
}

void LuisaDetector345::compute_denominator(int len) {
  derive_vector(&rho_to_derive_two, len, 2); // drho
  derive_vector(&rho_to_derive_two, len, 2); // d2rho
  derive_vector(&rho_to_derive_two, len, 2); // d3rho
  for (int i = 0; i < len; i++) {
    denominator[i] = fabs(rho_to_derive_two[i]);
  }
  derive_vector(&rho_to_derive_two, len, 2); // d4rho
  for (int i = 0; i < len; i++) {
    denominator[i] += fabs(rho_to_derive_two[i]);
  }
  derive_vector(&rho_to_derive_two, len, 2); // d5rho
  for (int i = 0; i < len; i++) {
    denominator[i] += fabs(rho_to_derive_two[i]);
  }
}

double LuisaDetector345::generic_first_der(int ind,
                                           const std::vector<double> &rho,
                                           int step) {
  auto l_flag = flag_functions::left(der_flags[ind]);
  auto r_flag = flag_functions::right(der_flags[ind]);
  if (std::min(l_flag, r_flag) >= 2 * step) {
    return 1 / 12. * (rho[ind - 2 * step] - rho[ind + 2 * step]) -
           2 / 3 * (rho[ind - step] - rho[ind + step]);
  }
  if (l_flag < 1 * step and r_flag >= 4 * step) {
    return -25 / 12. * rho[ind] + 4. * rho[ind + step] -
           3. * rho[ind + 2 * step] + 4 / 3. * rho[ind + 3 * step] -
           1 / 4. * rho[ind + 4 * step];
    /*return -11 / 6. * rho[ind] + 3 * rho[ind + step]
        - 3 / 2. * rho[ind + 2 * step] + 1 / 3. * rho[ind + 3 * step];*/
  }
  if (l_flag < 2 * step and r_flag >= 3 * step) {
    return -1 / 4. * rho[ind - step] - 5 / 6. * rho[ind] +
           3 / 2. * rho[ind + step] - 1 / 2. * rho[ind + 2 * step] +
           1 / 12. * rho[ind + 3 * step];
  }
  if (r_flag < 1 * step and l_flag >= 4 * step) {
    return 25 / 12. * rho[ind] - 4. * rho[ind - step] +
           3. * rho[ind - 2 * step] - 4 / 3. * rho[ind - 3 * step] +
           1 / 4. * rho[ind - 4 * step];
    /*return +11 / 6. * rho[ind] - 3 * rho[ind - step]
        + 3 / 2. * rho[ind - 2 * step] - 1 / 3. * rho[ind - 3 * step];*/
  }
  if (r_flag < 2 * step and l_flag >= 3 * step) {
    return +1 / 4. * rho[ind + step] + 5 / 6. * rho[ind] -
           3 / 2. * rho[ind - step] + 1 / 2. * rho[ind - 2 * step] -
           1 / 12. * rho[ind - 3 * step];
  }
  return 0.0;
}
