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

namespace PolySwab {
    // Source paper: An Online Algorithm for Segmenting Time Series
    // Source path: src/piecewise-approximation/polynomial/poly-swab.cpp

    class Approximator {
        private:
            static std::vector<long double> __interpolate(int degree, std::vector<long double>& data);
            static std::vector<long double> __regression(int degree, std::vector<long double>& data);
        
        public:
            static Polynomial approximate(std::string mode, int degree, std::vector<long double>& data);
            static long double cal_error(std::vector<long double>& segment, Polynomial& model);
    };

    class Grouper {
        public:
            static void merge(std::vector<long double>& s1, std::vector<long double>& s2, std::string mode);
            static long double merge_cost(std::vector<long double>& s1, std::vector<long double>& s2, std::string mode, int degree, long double error);
            static bool bound_check(std::vector<long double>& segment, Polynomial& model, long double error);
    };

    class Compression : public BaseCompression {
        private:
            int degree = -1;
            int n_segment = 0;
            long double error = 0;
            std::string mode = "";

            bool first = true;
            std::vector<long double> window;
            std::vector<std::vector<long double>> segments;

            void __bottom_up();
            bool __sliding_window();

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