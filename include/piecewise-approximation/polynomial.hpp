#include "base-c.hpp"
#include "base-d.hpp"


namespace CachedNormalEquation {
    // Source paper: Fast Piecewise Polynomial Fitting of Time-Series Data for Streaming Computing
    // Source path: src/piecewise-approximation/polynomial/normal-equation.cpp
    
    struct Model {
        int length = 0;
        float error = -1;
        Polynomial* function = nullptr;

        Model(Polynomial* function);
        ~Model();
    };

    class Compression : public BaseCompression {
        private:
            int degree = -1;
            long double error = 0;
            std::string mode = "";

            Model* model = nullptr;
            std::vector<Point2D> window;

            // std::map<int, Matrix<long double>*> cache;    // Use our matrix implementation
            std::map<int, Eigen::MatrixXd> cache;       // Use eigen library
            
            Polynomial* __calPolynomial();
            bool __approxSuccess(Model* model);

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            int degree = 0;

        protected:
            CSVObj* decompress(BinObj* compress_data) override;
        
        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace Swab {
    // Source paper: Fast Piecewise Polynomial Fitting of Time-Series Data for Streaming Computing
    // Source path: src/piecewise-approximation/polynomial/normal-equation.cpp
    struct Segment {
        std::vector<long double> window;
        std::vector<long double> coeffs;

        Segment(std::vector<long double> window, std::vector<long double> coeffs);
    };

    class Compression : public BaseCompression {
        private:
            int degree = -1;
            int n_segment = 0;
            long double error = 0;
            std::string mode = "";

            bool first = true;
            Segment* com_seg = nullptr;
            std::vector<long double> coeffs;
            std::vector<long double> window;
            std::vector<Segment> segments;

            std::vector<long double> __approximate(std::vector<long double> data);
            long double __verify(std::vector<long double>& segment, std::vector<long double>& coeffs);
            long double __merge_cost(Segment& s1, Segment& s2);
            void __merge(Segment& s1, Segment& s2);
            void __bottom_up();

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            int degree = -1;
            std::string mode = "";
            long double pivot = INFINITY;
            
        protected:
            CSVObj* decompress(BinObj* compress_data) override;
        
        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};