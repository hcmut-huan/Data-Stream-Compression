#include "base-c.hpp"
#include "base-d.hpp"

namespace Gorilla {
    // Source paper: Gorilla: A Fast, Scalable, In-Memory Time Series Database
    // Source path: src/floating-point/lossless/gorilla.cpp

    inline float END_SIGNAL = INFINITY;

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            BinObj* bitstream = nullptr;

            bool first = true;
            uint32_t storedVal;
            uint32_t storedLeadingZeros = INT_MAX;
            uint32_t storedTrailingZeros = INT_MAX;
            
            void __compress(uint32_t value);

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
            int buffer_size = 0;
            bool first = true;
            uint32_t storedVal;
            uint32_t storedLeadingZeros = INT_MAX;
            uint32_t storedTrailingZeros = INT_MAX;

            float __decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace Chimp {
    // Source paper: Chimp: Efficient Lossless Floating Point Compression for Time Series Databases
    // Source path: src/floating-point/lossless/chimp.cpp

    inline uint32_t END_SIGNAL = std::bit_cast<uint32_t>(INFINITY);

    class Compression : public BaseCompression {
        private:
            long size = 0;    
            int count = 0;
            int buffer_size = 0;
            std::string mode = "";
            BinObj* bitstream = nullptr;
            void (Compression::*__compress) (uint32_t);

            bool first = true;
            int setLsb;
            int* indices;
            int index = 0;
            int current = 0;
            int threshold;
            int previousValues;
            int previousValuesLog2;
            int storedLeadingZeros = INT_MAX;
            uint32_t storedVal;
            uint32_t* storedValues;
            
            short leadingRepresentation[32] = {
                0, 0, 0, 0, 0, 0, 1, 1,
                1, 1, 2, 2, 3, 3, 4, 4,
                5, 5, 6, 6, 7, 7, 7, 7,
                7, 7, 7, 7, 7, 7, 7, 7
            };

            short leadingRound[32] = {
                0, 0, 0, 0, 0, 0, 6, 6,
                6, 6, 10, 10, 12, 12, 14, 14,
                16, 16, 18, 18, 20, 20, 20, 20,
                20, 20, 20, 20, 20, 20, 20, 20
            };

            void __basic_compress(uint32_t value);
            void __128_compress(uint32_t value);

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
            int buffer_size = 0;
            std::string mode = "";
            float (Decompression::*__decompress) (BinObj*);

            bool first = true;
            int current = -1;
            int previousValues;
            int previousValuesLog2;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = 0;
            uint32_t storedVal;
            uint32_t* storedValues;
            
            short leadingRepresentation[8] = {0, 6, 10, 12, 14, 16, 18, 20};

            float __basic_decompress(BinObj* compress_data);
            float __128_decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace Elf {
    // Source paper: Elf: Erasing-based Lossless Floating-Point Compression
    // Source path: src/floating-point/lossless/elf.cpp

    inline uint32_t END_SIGNAL = std::bit_cast<uint32_t>(INFINITY);

    class Elf32Utils {
        private:
            static int f[21];
            static float map10iP[21];
            static float map10iN[21];
            static int mapSPGreater1[10];
            static float mapSPLess1[11];
            static double LOG_2_10;

        public:
            static int getFAlpha(int alpha);
            static std::pair<int, int> getAlphaAndBetaStar(float v, int lastBetaStar);
            static float roundUp(float v, int alpha);
            static int getSignificantCount(float v, int sp, int lastBetaStar);
            static float get10iP(int i);
            static float get10iN(int i);
            static int getSP(double v);
            static std::pair<int, int> getSPAnd10iNFlag(double v);
    };

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            BinObj* bitstream = nullptr;

            bool first = true;
            uint32_t storedVal = 0;
            int lastBetaStar = INT_MAX;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;

            short leadingRepresentation[32] = {
                0, 0, 0, 0, 0, 0, 1, 1,
                1, 1, 2, 2, 3, 3, 4, 4,
                5, 5, 6, 6, 7, 7, 7, 7,
                7, 7, 7, 7, 7, 7, 7, 7
            };

            short leadingRound[32] = {
                0, 0, 0, 0, 0, 0, 6, 6,
                6, 6, 10, 10, 12, 12, 14, 14,
                16, 16, 18, 18, 20, 20, 20, 20,
                20, 20, 20, 20, 20, 20, 20, 20
            };

            void __xorCompression(uint32_t value);
            void __compress(float value);

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
            bool first = true;
            int buffer_size = 0;
            uint32_t storedVal = 0;
            int lastBetaStar = INT_MAX;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;
            short leadingRepresentation[8] = {0, 6, 10, 12, 14, 16, 18, 20};

            float __recoverVByBetaStar(BinObj* compress_data);
            float __xorDecompress(BinObj* compress_data);
            float __decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace SElf {
    // Source paper: Adaptive Encoding Strategies for Lossless Floating-Point Compression
    // Source path: src/floating-point/lossless/self.cpp

    inline uint32_t END_SIGNAL = std::bit_cast<uint32_t>(INFINITY);

    struct Code {
        int code;
        int length;

        Code();
        Code(int code, int length);
    };

    struct Node {
        int data;
        int frequency;
        Node** children;
        int height; // leaf node is 0
        
        Node();
        ~Node();
        Node(int data);
        Node(int data, int frequency, int height);
    };

    struct CompareNode {
        bool operator()(const Node* a, const Node* b) const {
            if (a->frequency != b->frequency)
                return a->frequency > b->frequency;

            return a->height > b->height;
        }
    };

    class HuffmanEncode {
        private:
            static void freeTree(Node* root);
            static void generateHuffmanCodes(std::vector<Code>& huffmanCodes, Node* root, int code, int length);
            static bool nodeCompare(Node& a, Node& b);

        public:
            static std::vector<Code> getHuffmanCodes(int* frequencies);
            static Node* buildHuffmanTree(std::vector<Code> huffmanCodes);
            // static void writeHuffmanCodes(BinObj* bitstream, const std::vector<Code>& huffmanCodes);
            // static void readHuffmanCodes(BinObj* bitstream, std::vector<Code>& codes);
    };

    class SElf32Utils {
        private:
            static int f[21];
            static float map10iP[21];
            static float map10iN[21];
            static int mapSPGreater1[10];
            static float mapSPLess1[11];
            static double LOG_2_10;

        public:
            static int getFAlpha(int alpha);
            static std::pair<int, int> getAlphaAndBetaStar(float v, int lastBetaStar);
            static float roundUp(float v, int alpha);
            static int getSignificantCount(float v, int sp, int lastBetaStar);
            static float get10iP(int i);
            static float get10iN(int i);
            static int getSP(float v);
            static std::pair<int, int> getSPAnd10iNFlag(float v);
    };

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            BinObj* bitstream = nullptr;

            short leadingRepresentation[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 2, 2, 2, 2,
                    3, 3, 3, 3, 3, 3, 3, 3,
                    3, 3, 3, 3, 3, 3, 3, 3,
            };
            short leadingRound[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    8, 8, 8, 8, 12, 12, 12, 12,
                    16, 16, 16, 16, 16, 16, 16, 16,
                    16, 16, 16, 16, 16, 16, 16, 16
            };
            short trailingRepresentation[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1
            };
            short trailingRound[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    16, 16, 16, 16, 16, 16, 16, 16,
                    16, 16, 16, 16, 16, 16, 16, 16
            };

            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;
            
            bool first = true;
            uint32_t storedVal = 0;
            int leadingBitsPerValue = 2;
            int trailingBitsPerValue = 1;
                
            bool isFirstBlock = true;
            int lastBetaStar = INT_MAX;
            int frequency[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
            std::vector<Code> huffmanCode;

            void __xorCompression(uint32_t value);
            void __addValueFirst(float value);
            void __addValueHuffman(float value);
            void __compress(float value);

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
            int buffer_size = 0;
            int lastBetaStar = INT_MAX;
            bool isFirstBlock = true;
            Node* root;
            int frequency[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
            
            bool first = true;
            uint32_t storedVal = 0;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;
            short leadingRepresentation[8] = {0, 8, 12, 16};
            short trailingRepresentation[8] = {0, 16};
            int leadingBitsPerValue = 2;
            int trailingBitsPerValue = 1;

            float __recoverVByBetaStar(BinObj* compress_data);
            float __xorDecompress(BinObj* compress_data);
            float __decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};