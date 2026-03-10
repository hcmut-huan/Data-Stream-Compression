#include "piecewise-approximation/polynomial.hpp"

namespace PolySwab {
    Segment::Segment(std::vector<long double> window, std::vector<long double> coeffs) {
        this->window = window;
        this->coeffs = coeffs;
    }

    // Begin: compression
    std::vector<long double> Compression::__approximate(std::vector<long double> data) {
        std::vector<long double> coeffs;
        if (this->mode == "interpolate") {
            std::vector<int> indices;
            // Get evenly spaced indices between starting and ending data points
            for (int i = 0; i <= this->degree; i++) {
                int idx = (long)i * (data.size() - 1) / this->degree;
                indices.push_back(idx);
            }

            Eigen::MatrixXd A(this->degree+1, this->degree+1);
            Eigen::VectorXd b(this->degree+1);

            for (int i = 0; i < indices.size(); i++) {
                b(i) = data[indices[i]];
                for (int j = 0; j <= this->degree; j++) {
                    A(i, j) = std::pow(indices[i], j);
                }
            }

            Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);
            for (int i=0; i<=this->degree; i++) coeffs.push_back(c(i));
        }
        else if (this->mode == "regression") {
            Eigen::MatrixXd A(data.size(), this->degree + 1);
            Eigen::VectorXd b(data.size());

            for (int i = 0; i < data.size(); ++i) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                b(i) = data[i];
                for (int j = 0; j <= this->degree; ++j) {
                    A(i, j) = std::pow(i, j);
                }
            }
            
            Eigen::VectorXd c = (((A.transpose() * A).inverse()) * A.transpose()) * b;
            for (int i=0; i<=this->degree; i++) coeffs.push_back(c(i));
        }        

        return coeffs;
    }

    long double Compression::__verify(std::vector<long double>& segment, std::vector<long double>& coeffs) {
        long double max_error = 0;
        Polynomial poly(this->degree, coeffs);
        // Check that the approximated model satisfies the max error constraint
        for (int i=0; i<segment.size(); i++) {
            long double error = std::abs(segment[i]-poly.subs(i));
            max_error = error > max_error ? error : max_error;
            
            // Immediately terminate when individual error exceed the allowable threshold
            if (max_error > this->error) INFINITY;
        }

        return max_error;
    }

    long double Compression::__merge_cost(Segment& s1, Segment& s2) {
        int s1_len = s1.window.size();
        std::vector<long double> s = s1.window;

        auto start = this->mode == "interpolate" ? s2.window.begin() + 1 : s2.window.begin();
        s.insert(s.end(), start, s2.window.end());
        // merge cost is related to the maximum individual error
        std::vector<long double> coeffs = this->__approximate(s);
        long double error = this->__verify(s, coeffs);

        return error;
    }

    void Compression::__merge(Segment& s1, Segment& s2) {
        int s1_len = s1.window.size();
        auto start = this->mode == "interpolate" ? s2.window.begin() + 1 : s2.window.begin();
        s1.window.insert(s1.window.end(), start, s2.window.end());
        s1.coeffs = this->__approximate(s1.window);
    }

    void Compression::__bottom_up() {
        std::vector<long double> m_err;
        for (int i=0; i<this->n_segment-1; i++) {
            m_err.push_back(this->__merge_cost(this->segments[i], this->segments[i+1]));
        }

        auto it = std::min_element(m_err.begin(), m_err.end());
        while (*it <= this->error && this->segments.size() > 1) {
            int index = std::distance(m_err.begin(), it);
            this->__merge(this->segments[index], this->segments[index+1]);

            m_err.erase(m_err.begin() + index);
            this->segments.erase(this->segments.begin() + index + 1);
            
            if (index != 0) m_err[index-1] = this->__merge_cost(this->segments[index-1], this->segments[index]);
            if (index < m_err.size()) m_err[index] = this->__merge_cost(this->segments[index], this->segments[index+1]);
            it = std::min_element(m_err.begin(), m_err.end());  
        }
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
        this->degree = atoi(params[2]);
        this->n_segment = atoi(params[3]);
    }

    void Compression::finalize() {
        if (this->window.size() >= this->degree + 1) {
            this->segments.push_back(Segment(this->window, this->coeffs));
        }

        if (this->segments.size() > 1) {
            this->__bottom_up();
        }

        for (int i=0; i<this->segments.size(); i++) {
            this->com_seg = &this->segments[i];
            this->yield();
        }
        this->segments.clear();
    }

    void Compression::compress(Univariate* data) {
        this->window.push_back(data->get_value());

        if (this->window.size() == this->degree + 1) {
            this->coeffs = this->__approximate(this->window);
        }
        else if (this->window.size() > this->degree + 1) {
            std::vector<long double> n_coeffs = this->__approximate(this->window);
            if (this->__verify(this->window, n_coeffs) <= this->error) {
                this->coeffs = n_coeffs;
            }
            else {
                this->window.pop_back();
                this->segments.push_back(Segment(this->window, this->coeffs));
                
                if (this->mode == "interpolate") {
                    // Add last data point of previous segment to ensure connectivity
                    long double tail = this->window.back();
                    this->window.clear();
                    this->window.push_back(tail);
                    this->window.push_back(data->get_value());

                    if (this->window.size() == this->degree + 1) {
                        this->coeffs = this->__approximate(this->window);
                    }
                }
                else if (this->mode == "regression") {
                    this->window.clear();
                    this->window.push_back(data->get_value());
                }
            }
            // Perform bottom up once n_segment sliding window segments have accumulated
            if (this->segments.size() == this->n_segment) {
                this->__bottom_up();
                this->com_seg = &this->segments[0];
                this->yield();
                this->com_seg = nullptr;
                this->segments.erase(this->segments.begin());
            }
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        if (this->mode == "interpolate") {
            if (this->first) {
                obj->put((float) this->com_seg->coeffs[0]);
                this->first = false;
            }
            
            obj->put(VariableByteEncoding::encode(this->com_seg->window.size()));
            for (int i = 1; i <= this->degree; i++) {
                obj->put((float) this->com_seg->coeffs[i]);
            }
        }
        else if (this->mode == "regression") {
            obj->put(VariableByteEncoding::encode(this->com_seg->window.size()));
            for (int i = 0; i <= this->degree; i++) {
                obj->put((float) this->com_seg->coeffs[i]);
            }
        }
        
        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        this->degree = atoi(params[2]);   
        this->mode = params[1];
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        if (this->mode == "interpolate") {
            int start = 1;
            if (this->pivot == INFINITY) {
                start = 0;
                this->pivot = compress_data->getFloat();
            }

            unsigned long length = VariableByteEncoding::decode(compress_data);
            float* coefficients = new float[this->degree+1];
            for (int i = 1; i <= this->degree; i++) {
                coefficients[i] = compress_data->getFloat();
            }
            coefficients[0] = this->pivot;
            Polynomial model(this->degree, coefficients);

            for (int i=start; i<length; i++) {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                    base_obj->pushData(std::to_string(model.subs(i)));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime + i * this->interval));
                    obj->pushData(std::to_string(model.subs(i)));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
            }

            delete[] coefficients;
            this->basetime += (length - 1) * this->interval;
            this->pivot = model.subs(length-1);
        }
        else if (this->mode == "regression") {
            unsigned long length = VariableByteEncoding::decode(compress_data);
            float* coefficients = new float[this->degree+1];
            for (int i = 0; i <= this->degree; i++) {
                coefficients[i] = compress_data->getFloat();
            }
            Polynomial model(this->degree, coefficients);

            for (int i=0; i<length; i++) {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                    base_obj->pushData(std::to_string(model.subs(i)));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime + i * this->interval));
                    obj->pushData(std::to_string(model.subs(i)));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
            }

            delete[] coefficients;
            this->basetime += length * this->interval;
        }

        
        return base_obj;
    }
    // End: decompression
};