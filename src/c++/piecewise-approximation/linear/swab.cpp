#include "piecewise-approximation/linear.hpp"

namespace Swab {
    Segment::Segment(std::vector<long double>& data, UpperHull& l_cvx, LowerHull& u_cvx) {
        this->data = data;
        this->u_cvx = u_cvx;
        this->l_cvx = l_cvx;
    }

    Line Approximator::__interpolate(std::vector<long double>& data) {
        Point2D first(0, data.front());
        Point2D last(data.size()-1, data.back());

        return Line::line(first, last);
    }

    Line Approximator::__regression(std::vector<long double>& data) {
        Eigen::MatrixXd A(data.size(), 2);
        Eigen::VectorXd b(data.size());
        for (int i = 0; i < data.size(); ++i) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
            b(i) = data[i];
            A(i, 0) = 1; A(i, 1) = i;
        }
        Eigen::VectorXd c = (((A.transpose() * A).inverse()) * A.transpose()) * b;
        
        return Line(c(1), c(0));
    }

    Line Approximator::approximate(std::string mode, std::vector<long double>& data) {
        if (mode == "regression")  return Approximator::__regression(data);
        else return Approximator::__interpolate(data);
    }

    // Mean square error calculate
    long double Approximator::cal_error(std::vector<long double>& data, Line& line) {
        long double error = 0;
        for (int i=0; i<data.size(); i++) {
            error += std::pow(data[i]-line.subs(i), 2);
        }

        return error / data.size();
    }

    // Verify new line satisfies infinity bound or not
    bool Grouper::bound_check(Segment& segment, Line& line, int offset) {
        UpperHull& l_cvx = segment.l_cvx;
        LowerHull& u_cvx = segment.u_cvx;
        
        // Immediately terminate when individual error exceed the allowable threshold
        for (int i=0; i<l_cvx.size(); i++) {
            Point2D p = l_cvx.at(i);
            if (line.subs(p.x + offset) - p.y < 0) return false;
        }

        for (int i=0; i<u_cvx.size(); i++) {
            Point2D p = u_cvx.at(i);
            if (p.y - line.subs(p.x + offset) < 0) return false;
        }

        return true;
    }

    void Grouper::merge(Segment& s1, Segment& s2, std::string mode) {
        int offset = mode == "interpolate" ? 1 : 0;
        s1.u_cvx.concat(s2.u_cvx, s1.data.size() - offset);
        s1.l_cvx.concat(s2.l_cvx, s1.data.size() - offset);
        s1.data.insert(s1.data.end(), s2.data.begin() + offset, s2.data.end());
    }

    long double Grouper::merge_cost(Segment& s1, Segment& s2, std::string mode) {
        int offset = mode == "interpolate" ? 1 : 0;
        std::vector<long double> s = s1.data;
        s.insert(s.end(), s2.data.begin() + offset, s2.data.end());
        Line line = Approximator::approximate(mode, s);
    
        if (!Grouper::bound_check(s1, line)) return INFINITY;
        else if (!Grouper::bound_check(s2, line, s1.data.size()-offset)) return INFINITY;
        else return Approximator::cal_error(s, line);
    }

    // Begin: compression
    bool Compression::__sliding_window() {
        Point2D p(this->window.size()-1, this->window.back());
        Line line = Approximator::approximate(this->mode, this->window);
        Segment segment(this->window, this->l_cvx, this->u_cvx);

        return Grouper::bound_check(segment, line) && 
            std::abs(line.subs(p.x) - p.y) <= this->error;
    }

    void Compression::__bottom_up() {
        std::vector<long double> m_err;
        for (int i=0; i<this->segments.size()-1; i++) {
            m_err.push_back(Grouper::merge_cost(
                this->segments[i], this->segments[i+1], this->mode
            ));
        }

        auto it = std::min_element(m_err.begin(), m_err.end());
        while (*it != INFINITY && this->segments.size() > 1) {
            int index = std::distance(m_err.begin(), it);
            Grouper::merge(this->segments[index], this->segments[index+1], this->mode);

            m_err.erase(m_err.begin() + index);
            this->segments.erase(this->segments.begin() + index + 1);
            
            if (index != 0) {
                m_err[index-1] = Grouper::merge_cost(
                    this->segments[index-1], this->segments[index], this->mode
                );
            }
            if (index < m_err.size()) {
                m_err[index] = Grouper::merge_cost(
                    this->segments[index], this->segments[index+1], this->mode
                );
            }
            
            it = std::min_element(m_err.begin(), m_err.end());  
        }
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
        this->n_segment = atoi(params[2]);
    }

    void Compression::finalize() {
        if (this->window.size() >= 2) {
            this->segments.push_back(Segment(
                this->window, this->l_cvx, this->u_cvx
            ));
        }

        this->__bottom_up();
        while (this->segments.size() != 0) {
            this->yield();
            this->segments.erase(this->segments.begin());
        }
    }

    void Compression::compress(Univariate* data) {
        this->window.push_back(data->get_value());

        if (this->window.size() > 2) {
            if (!this->__sliding_window()) {
                this->window.pop_back();
                this->segments.push_back(Segment(
                    this->window, this->l_cvx, this->u_cvx
                ));
                
                if (this->mode == "interpolate") {
                    // Add last data point of previous segment to ensure connectivity
                    long double tail = this->window.back();
                    this->window.clear();
                    this->window.push_back(tail);
                    this->window.push_back(data->get_value());

                    this->l_cvx.clear();
                    this->l_cvx.append(Point2D(0, tail - this->error));
                    this->u_cvx.clear();
                    this->u_cvx.append(Point2D(0, tail + this->error));
                }
                else if (this->mode == "regression") {
                    this->window.clear();
                    this->window.push_back(data->get_value());
                    this->l_cvx.clear();
                    this->u_cvx.clear();
                }
            }
        }

        // Perform bottom up once n_segment sliding window segments have accumulated
        if (this->segments.size() == this->n_segment) {
            this->__bottom_up();
            this->yield();
            this->segments.erase(this->segments.begin());
        }

        this->l_cvx.append(Point2D(this->window.size()-1, data->get_value() - this->error));
        this->u_cvx.append(Point2D(this->window.size()-1, data->get_value() + this->error));
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        Segment& segment = this->segments[0];

        if (this->mode == "interpolate") {
            if (this->first) {
                obj->put((float) segment.data.front());
                this->first = false;
            }
            
            obj->put(VariableByteEncoding::encode(segment.data.size()));
            obj->put((float) segment.data.back());
        }
        else if (this->mode == "regression") {
            Line line = Approximator::approximate(this->mode, segment.data);
            obj->put(VariableByteEncoding::encode(segment.data.size()));
            obj->put((float) line.get_intercept());
            obj->put((float) line.get_slope());
        }
        
        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        this->mode = params[1];
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int start = 0;
        unsigned long length = 0;
        Line line(0, 0);
        int offset = this->mode == "interpolate" ? 1 : 0;

        if (this->mode == "interpolate") {
            start = 1;
            if (this->pivot == INFINITY) {
                start = 0;
                this->pivot = compress_data->getFloat();
            }

            length = VariableByteEncoding::decode(compress_data);
            float last = compress_data->getFloat();
            line = Line::line(Point2D(0, this->pivot), Point2D(length-1, last)); 
            this->pivot = last;
        }
        else if (this->mode == "regression") {
            length = VariableByteEncoding::decode(compress_data);
            float intercept = compress_data->getFloat();
            float slope = compress_data->getFloat();
            line = Line(slope, intercept);
        }

        for (int i=start; i<length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(line.subs(i)));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(line.subs(i)));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }
        
        this->basetime += (length - offset) * this->interval;
        return base_obj;
    }
    // End: decompression
};