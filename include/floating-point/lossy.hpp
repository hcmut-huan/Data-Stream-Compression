#include "base-c.hpp"
#include "base-d.hpp"

// namespace Camel {
//     // Source paper: Camel: Efficient Compression of Floating-Point Time Series
//     // Source path: src/floating-point/lossy/camel.cpp

//     class Compression : public BaseCompression {
//         private:
//             long buffer = 0;

//         protected:
//             void compress(Univariate* data) override;
//             BinObj* serialize() override;

//         public:
//             Compression(std::string output) : BaseCompression(output) {}            
//             void initialize(int count, char** params) override;
//             void finalize() override;
//     };

//     class Decompression : public BaseDecompression {
//         protected:
//             CSVObj* decompress(BinObj* compress_data) override;

//         public:
//             Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
//             void initialize(int count, char** params) override;
//             void finalize() override;
//     };
// };

namespace Serf {
    // Source paper: Serf: Streaming Error-Bounded Floating-Point Compression
    // Source path: src/floating-point/lossy/serf.cpp

    inline float END_SIGNAL = INFINITY;

    class SerfUtils32 {
        public:
            static uint32_t FindAppInt(float min, float max, float v, uint32_t last_int, float max_diff);

        private:
            static constexpr uint32_t kBitWeight[32] = {
                1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U, 256U, 512U, 1024U,
                2048U, 4096U, 8192U, 16384U, 32768U, 65536U, 131072U, 262144U,
                524288U, 1048576U, 2097152U, 4194304U, 8388608U, 16777216U, 33554432U,
                67108864U, 134217728U, 268435456U, 536870912U, 1073741824U, 2147483648U
            };

            static uint32_t FindAppInt(float min_float, float max_float, uint32_t sign, 
                float original, uint32_t last_int, float max_diff);
    };

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            long double error;
            int buffer_size = 0;
            std::string mode = "";
            BinObj* bitstream = nullptr;
            void (Compression::*__compress) (float);

            bool first = true;
            float storedVal = 2;
            int leading_bits_per_value_ = 2;
            int trailing_bits_per_value_ = 1;
            int lead_distribution_[32];
            int trail_distribution_[32];
            int stored_leading_zeros_ = INT_MAX;
            int stored_trailing_zeros_ = INT_MAX;

            long compressed_size_this_block_;
            long compressed_size_last_block_ = 0;
            long compressed_size_this_window_ = 0;
            int number_of_values_this_window_ = 0;
            double compression_ratio_last_window_ = 0;


            short leading_representation_[32] = {
                0, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 2, 2, 2, 2,
                3, 3, 3, 3, 3, 3, 3, 3,
                3, 3, 3, 3, 3, 3, 3, 3,
            };

            short leading_round_[32] = {
                0, 0, 0, 0, 0, 0, 0, 0,
                8, 8, 8, 8, 12, 12, 12, 12,
                16, 16, 16, 16, 16, 16, 16, 16,
                16, 16, 16, 16, 16, 16, 16, 16
            };

            short trailing_representation_[32] = {
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1
            };

            short trailing_round_[32] = {
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                16, 16, 16, 16, 16, 16, 16, 16,
                16, 16, 16, 16, 16, 16, 16, 16
            };

            void __qt_compress(float value);
            void __xor_compress(float value);

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
            
            long double error;
            int buffer_size = 0;
            std::string mode = "";
            float (Decompression::*__decompress) (BinObj*);

            bool first = true;
            float storedVal = 2;
            int stored_leading_zeros_ = INT_MAX;
            int stored_trailing_zeros_ = INT_MAX;
            short leading_representation_[4] = {0, 8, 12, 16};
            short trailing_representation_[2] = {0, 16};
            int leading_bits_per_value_ = 2;
            int trailing_bits_per_value_ = 1;
            short look_up_table[16] = {16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

            float __qt_decompress(BinObj* compress_data);
            float __xor_decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};