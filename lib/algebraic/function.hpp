#ifndef ALGEBRAIC_FUNCTION_HPP
#define ALGEBRAIC_FUNCTION_HPP

#include <string>
#include <cmath>

struct Point2D {
    long double x;
    long double y;

    Point2D() {
        this->x = 0;
        this->y = 0;
    }

    Point2D(long double x, long double y) {
        this->x = x;
        this->y = y;
    }
};

// Geometry function support timeseries analysis
class Line {
    private:
        long double slope;
        long double intercept;

    public:
        Line() {
            this->slope = 0;
            this->intercept = 0;
        }

        Line(long double slope, long double intercept) {
            this->slope = slope;
            this->intercept = intercept;
        }

        long double subs(long double x) {
            return this->slope * x + this->intercept;
        }

        void set_slope(long double slope) {
            this->slope = slope;
        }

        long double get_slope() {
            return this->slope;
        }

        long double get_intercept() {
            return this->intercept;
        }

        void set_intercept(long double intercept) {
            this->intercept = intercept;
        }

        long double get_root() {
            if (this->slope == 0) return INFINITY;
            else return - this->intercept / this->slope;
        }

        static Line line(Point2D p1, Point2D p2) {
            long double slope = (p1.y - p2.y) / (p1.x - p2.x);
            long double intercept = p1.y - slope*p1.x;

            return Line(slope, intercept);
        }

        static Line line(long double slope, Point2D p) {
            return Line(slope, p.y - slope*p.x);
        }

        static Point2D intersection(Line l1, Line l2) {
            long double x = (l1.intercept - l2.intercept) / (l2.slope - l1.slope);
            long double y = l1.intercept + l1.slope*x;

            return Point2D(x, y); 
        }

};

// Source for manipulating polynomial function
class Polynomial {
    public:
        int degree;
        long double* coefficients;    // coefficient degree starts from 0 

    public:
        Polynomial(int k, const float* coeffs, bool reverse = false) {
            this->degree = k;
            this->coefficients = new long double[k+1];
            if (reverse) for (int i=0; i<k+1; i++) this->coefficients[k-i] = coeffs[i];
            else for (int i=0; i<k+1; i++) this->coefficients[i] = coeffs[i];
        }

        Polynomial(int k, const double* coeffs, bool reverse = false) {
            this->degree = k;
            this->coefficients = new long double[k+1];
            if (reverse) for (int i=0; i<k+1; i++) this->coefficients[k-i] = coeffs[i];
            else for (int i=0; i<k+1; i++) this->coefficients[i] = coeffs[i];
        }

        Polynomial(int k, const long double* coeffs, bool reverse = false) {
            this->degree = k;
            this->coefficients = new long double[k+1];
            if (reverse) for (int i=0; i<k+1; i++) this->coefficients[k-i] = coeffs[i];
            else for (int i=0; i<k+1; i++) this->coefficients[i] = coeffs[i];
        }

        Polynomial(float coeff) {
            this->degree = 0;
            this->coefficients = new long double[1];
            this->coefficients[0] = (long double) coeff;
        }

        Polynomial(long double coeff) {
            this->degree = 0;
            this->coefficients = new long double[1];
            this->coefficients[0] = coeff;
        }

        ~Polynomial() {
            delete[] this->coefficients;
        }

        long double subs(long double indeterminate) const {
            long double result = this->coefficients[0];
            for (int i=1; i<this->degree+1; i++) {
                result += this->coefficients[i]*pow(indeterminate, i);
            }

            return result;
        }

        long double get_coefficient(int degree) {
            return this->coefficients[degree];
        }

        std::string str() const {
            std::string s = "";
            for (int i=this->degree; i>0; i--) {
                s += std::to_string(this->coefficients[i]) + "x^" + std::to_string(i) + " + ";
            }

            return s + std::to_string(this->coefficients[0]);
        }
};

#endif