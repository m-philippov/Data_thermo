#include "effective.h"


Effective::Effective(const std::vector<float>& T, const float& Tcold, const float& Thot, const std::vector<float>& X, 
    const std::vector<float>& Y, const float& characteristic_size, const std::string& Xname, const std::string& Yname, const std::string& filename)
{
    for (int i = 0; i < filename.size() - 4; ++i) {
        _filename += filename[i];
    }
    _filename += "eff.txt";
    _characteristic_size = characteristic_size;
    _X.resize(X.size()); 
    for (int i = 0; i < X.size(); i++) {
        _X[i] = X[i]/characteristic_size;
    }
    _Y.resize(Y.size());
    for (int i = 0; i < Y.size(); i++) {
        _T[i] = Y[i]/characteristic_size;
    }
   
    _X_Name = Xname;
    _Y_Name = Yname;
    _T.resize(T.size());
    for (int i = 0; i < T.size(); i++) {
        _T[i] = T[i] ;
    }
    _eff.resize(T.size());
#pragma omp parallel for
    for (int i = 0; i < T.size(); i++) {
        _eff[i] = (T[i] - Tcold) / (Thot - Tcold);
    }
}

void Effective::PrintEffective(const std::string& delimetr)
{
    std::ofstream f;
    f.open(_filename, std::ios::out);
    f << _X_Name << delimetr << _Y_Name << delimetr << "T" << delimetr << "eff" << '\n';

    for (int i = 0; i < _T.size(); i++) {
        f << _X[i] << delimetr << _Y[i] << delimetr << _T[i] << delimetr << _eff[i] << '\n';
    }
    f.close();
}

void AreaEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name)
{
    float Ti;
    for (int i = 0; i < bC.size(); i++) {
        std::vector<float> T_;
        std::vector<float> X;
        std::vector<float> Y;
        std::vector<float> eff;
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            infile >> Ti;
            if ((j % width_pix + margin_from_edges_pix < bC[i].x1) && (j % width_pix - margin_from_edges_pix > bC[i].x0) && (j / width_pix - margin_from_edges_pix > bC[i].ymax) && (j / width_pix + margin_from_edges_pix < bC[i].ymin)) {
                T_.push_back(linear_ax * Ti - linear_b);
                X.push_back(((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
                Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
                eff.push_back((bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1));
            }
        }
        infile.close();
        std::string _filename = "eff_area_" + bC[i].number + ".txt";

        std::ofstream f;
        f.open(_filename, std::ios::out);
        f << x_name << ", " << y_name << ", " << temp_name << ", " << eff_name << '\n';
#pragma omp parallel for
        for (int i = 0; i < T_.size(); i++) {
            f << X[i] << ", " << Y[i] << ", " << T_[i] << ", " << eff[i] << '\n';
        }
        f.close();
    }
}

void AreaIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name) {
    float Ti;
    std::vector<float> T_;
    std::vector<float> X;
    std::vector<float> Y;
    std::vector<float> eff;
    for (int i = 0; i < bC.size(); i++) {
        float k = static_cast<float>(0.0001);
        T_.push_back(0);
        X.push_back(0);
        Y.push_back(0);
        eff.push_back(0);
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            infile >> Ti;
            if ((j % width_pix + margin_from_edges_pix < bC[i].x1) && (j % width_pix - margin_from_edges_pix > bC[i].x0) && (j / width_pix - margin_from_edges_pix > bC[i].ymax) && (j / width_pix + margin_from_edges_pix < bC[i].ymin)) {
                k++;
                T_[i] += linear_ax * Ti - linear_b;
                X[i] += ((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff;
                Y[i] += (bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff;
                eff[i] += (bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1);
            }
        }
        T_[i] /= k;
        X[i] /= k;
        Y[i] /= k;
        eff[i] /= k;
        infile.close();
       
    }
    std::ofstream f;
    std::string _filename = "eff_area_integral.txt";
    f.open(_filename, std::ios::out);
    f << "num" << "\t" << x_name << "\t" << y_name << "\t" << temp_name << "\t" << eff_name << '\n';
#pragma omp parallel for
    for (int i = 0; i < T_.size(); i++) {
        f << bC[i].number << "\t" << X[i] << "\t" << Y[i] << "\t" << T_[i] << "\t" << eff[i] << '\n';
    }
    f.close();
}
void AlongXbEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name) {
    for (int i = 0; i < bC.size(); i++) {
        float Ti;
        std::vector<float> T_;
        std::vector<float> X;
        std::vector<float> Y;
        std::vector<float> eff;
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            infile >> Ti;
            if ((j % width_pix == bC[i].xb) && (j / width_pix - margin_from_edges_pix > bC[i].ymax) && (j / width_pix + margin_from_edges_pix < bC[i].ymin)) {
                T_.push_back(linear_ax * Ti - linear_b);
                X.push_back(((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
                Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
                eff.push_back((bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1));
            }
        }
        infile.close();
        std::string _filename = "eff_alongX_" + bC[i].number + ".txt";

        std::ofstream f;
        f.open(_filename, std::ios::out);
        f << y_name << "\t" << temp_name << "\t" << eff_name << '\n';
#pragma omp parallel for
        for (int i = 0; i < T_.size(); i++) {
            f  << Y[i] << "\t" << T_[i] << "\t" << eff[i] << '\n';
        }
        f.close();
    }
}

void AlongXEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float X0) {
    for (int i = 0; i < bC.size(); i++) {
        float Ti;
        std::vector<float> T_;
        std::vector<float> X;
        std::vector<float> Y;
        std::vector<float> eff;
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            infile >> Ti;
            if ((j % width_pix == static_cast<int>(bC[i].xb + (X0 + x_begin_coord)/ bC[i].dx_mm)) && (j / width_pix - margin_from_edges_pix > bC[i].ymax) && (j / width_pix + margin_from_edges_pix < bC[i].ymin)) {
                T_.push_back(linear_ax * Ti - linear_b);
                X.push_back(((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
                Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
                eff.push_back((bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1));
            }
        }
        infile.close();
        std::string _filename = "eff_alongX_" + std::to_string(static_cast<int>(X0)) + "mm" + bC[i].number  + ".txt";
        std::ofstream f;
        f.open(_filename, std::ios::out);
        f << y_name << "\t" << temp_name << "\t" << eff_name << '\n';
#pragma omp parallel for
        for (int i = 0; i < T_.size(); i++) {
            f << Y[i] << "\t" << T_[i] << "\t" << eff[i] << '\n';
        }
        f.close();
    }
}
void AlongYEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float y0) {
    for (int i = 0; i < bC.size(); i++) {
        float Ti;
        std::vector<float> T_;
        std::vector<float> X;
        std::vector<float> Y;
        std::vector<float> eff;
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            infile >> Ti;
            if ((j % width_pix + margin_from_edges_pix < bC[i].x1) && (j % width_pix - margin_from_edges_pix > bC[i].x0) && (j / width_pix - margin_from_edges_pix == static_cast<int>((y0 + y_begin_coord) / bC[i].dx_mm))) {
                T_.push_back(linear_ax * Ti - linear_b);
                X.push_back(((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
                Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
                eff.push_back((bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1));
            }
        }
        infile.close();
        std::string _filename = "eff_alongY_" + std::to_string(static_cast<int>(y0)) + "mm" + bC[i].number + ".txt";
        std::ofstream f;
        f.open(_filename, std::ios::out);
        f << x_name << "\t" << temp_name << "\t" << eff_name << '\n';
#pragma omp parallel for
        for (int i = 0; i < T_.size(); i++) {
            f << X[i] << "\t" << T_[i] << "\t" << eff[i] << '\n';
        }
        f.close();
    }
}

void AlongXbIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float dx) {
    float Ti;
    for (int i = 0; i < bC.size(); i++) {
        std::vector<float> T_;
        std::vector<float> X;
        std::vector<float> Y;
        std::vector<float> eff;
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            int i = 0;
            int k = 0;
            float p = FLT_MIN;
            infile >> Ti;
            if ((j % width_pix + margin_from_edges_pix < bC[i].xb + dx / bC[i].dx_mm/2) && (j % width_pix - margin_from_edges_pix > bC[i].xb - dx / bC[i].dx_mm / 2) && (j / width_pix - margin_from_edges_pix > bC[i].ymax) && (j / width_pix + margin_from_edges_pix < bC[i].ymin)) {
                
                if ((p <= j / width_pix - 0.0000001) || (p >= j / width_pix + 0.0000001)) {
                    if (i != 0) {
                        T_[i] /= k;
                        eff[i] /= k;
                    }
                    p = (((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
                    i++;
                    T_.push_back(linear_ax * Ti - linear_b);
                    Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
                    eff.push_back((bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1));
                    k = 1;
                }
                else {
                    T_[i] += linear_ax * Ti - linear_b;
                    eff[i] += (bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1);
                    k++;
                }
            }
        }
        infile.close();
        std::string _filename = "eff_area_integral_Xb" + std::to_string(static_cast<int>(dx)) + bC[i].number + ".txt";

        std::ofstream f;
        f.open(_filename, std::ios::out);
        f << y_name << "\t" << temp_name << "\t" << eff_name << '\n';
#pragma omp parallel for
        for (int i = 0; i < T_.size(); i++) {
            f << Y[i] << "\t" << T_[i] << "\t" << eff[i] << '\n';
        }
        f.close();
    }
}
void AlongXIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float dx, float x0) {
    float Ti;
    for (int i = 0; i < bC.size(); i++) {
        std::vector<float> T_;
        std::vector<float> X;
        std::vector<float> Y;
        std::vector<float> eff;
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            int i = 0;
            int k = 0;
            float p = FLT_MIN;
            infile >> Ti;
            if ((j % width_pix + margin_from_edges_pix < bC[i].xb + (x0 + x_begin_coord + dx / 2) / bC[i].dx_mm) && (j % width_pix - margin_from_edges_pix > bC[i].xb + (x0 + x_begin_coord - dx / 2) / bC[i].dx_mm) && (j / width_pix - margin_from_edges_pix > bC[i].ymax) && (j / width_pix + margin_from_edges_pix < bC[i].ymin)) {

                if ((p <= (((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff) - 0.00001) || (p >= (((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff) + 0.00001)) {
                    p = (((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
                    if (i != 0) {
                        T_[i] /= k;
                        eff[i] /= k;
                    }
                    i++;
                    T_.push_back(linear_ax * Ti - linear_b);
                    Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
                    eff.push_back((bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1));
                    k = 1;
                }
                else {
                    T_[i] += linear_ax * Ti - linear_b;
                    eff[i] += (bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1);
                    k++;
                }
            }
        }
        infile.close();
        std::string _filename = "eff_area_integral_X" + std::to_string(static_cast<int>(dx)) + bC[i].number + ".txt";

        std::ofstream f;
        f.open(_filename, std::ios::out);
        f << y_name << "\t" << temp_name << "\t" << eff_name << '\n';
#pragma omp parallel for
        for (int i = 0; i < T_.size(); i++) {
            f << Y[i] << "\t" << T_[i] << "\t" << eff[i] << '\n';
        }
        f.close();
    }
}
void AlongYIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float dy, float y0) {
    float Ti;
    for (int i = 0; i < bC.size(); i++) {
        std::vector<float> T_;
        std::vector<float> X;
        std::vector<float> Y;
        std::vector<float> eff;
        std::ifstream infile(bC[i].number + ".txt");
#pragma omp parallel for
        for (int j = 0; j < height_pix * width_pix; j++) {
            int i = 0;
            int k = 0;
            float p = FLT_MIN;
            infile >> Ti;
            if ((j % width_pix + margin_from_edges_pix < bC[i].x1) && (j % width_pix - margin_from_edges_pix > bC[i].x0) && (j / width_pix - margin_from_edges_pix > (y0 + y_begin_coord + dy / 2) / bC[i].dx_mm) && (j / width_pix + margin_from_edges_pix < (y0 + y_begin_coord - dy / 2) / bC[i].dx_mm)) {

                if ((p <= (((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff) - 0.00001) || (p >= (((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff) + 0.00001)) {
                    p = (((j % width_pix) * bC[i].dx_mm - bC[i].xb * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
                    if (i != 0) {
                        T_[i] /= k;
                        eff[i] /= k;
                    }
                    i++;
                    T_.push_back(linear_ax * Ti - linear_b);
                    Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
                    eff.push_back((bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1));
                    k = 1;
                }
                else {
                    T_[i] += linear_ax * Ti - linear_b;
                    eff[i] += (bC[i].t3 - (linear_ax * Ti - linear_b)) / (bC[i].t3 - bC[i].t1);
                    k++;
                }
            }
        }
        infile.close();
        std::string _filename = "eff_area_integral_Y" + std::to_string(static_cast<int>(dy)) + bC[i].number + ".txt";

        std::ofstream f;
        f.open(_filename, std::ios::out);
        f << x_name << "\t" << temp_name << "\t" << eff_name << '\n';
#pragma omp parallel for
        for (int i = 0; i < T_.size(); i++) {
            f << X[i] << "\t" << T_[i] << "\t" << eff[i] << '\n';
        }
        f.close();
    }
}
