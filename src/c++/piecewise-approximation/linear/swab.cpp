#include "piecewise-approximation/linear.hpp"

namespace Swab {
    // Begin: Base segment implementation
    long double Swab::__verify(
        UpperHull& l_cvx, LowerHull& u_cvx, 
        Line& n_line, long double threshold, int offset) 
    {
        long double max_error = 0;
        // Check that the approximated model satisfies the max error constraint
        for (int i=0; i<l_cvx.size(); i++) {
            Point2D p = l_cvx.at(i);
            long double error = n_line.subs(p.x+offset) - p.y;
            
            // Immediately terminate when individual error exceed the allowable threshold
            if (error < 0) return INFINITY;
            else max_error = error > max_error ? error : max_error;
        }

        for (int i=0; i<u_cvx.size(); i++) {
            Point2D p = u_cvx.at(i);
            long double error = p.y-n_line.subs(p.x+offset);

            if (error < 0) return INFINITY;
            else max_error = error > max_error ? error : max_error;
        }

        return std::abs(threshold - max_error);
    }
    // End: base segment implementation

    // Begin: interpolate segment implementation
    Line InterpolateSegment::__approximate(Point2D* first, Point2D* last, int offset) {
        Point2D n_last(last->x + offset, last->y);
        return Line::line(*first, n_last);
    }

    InterpolateSegment::~InterpolateSegment() {
        this->u_cvx.clear();
        this->l_cvx.clear();
        if (this->first != nullptr) delete this->first;
        if (this->last != nullptr) delete this->last;
        if (this->line != nullptr) delete this->line;
    }

    BinObj* InterpolateSegment::serialize(bool& first) {
        BinObj* obj = new BinObj;
        if (first) {
            obj->put((float) this->first->y);
            first = false;
        }
        
        obj->put(VariableByteEncoding::encode(this->size));
        obj->put((float) this->last->y);

        // std::cout << "compress: " << this->size << " " << this->line->get_slope() << " " << this->line->get_intercept() << "\n";

        return obj;
    }

    long double InterpolateSegment::approximate(long double data, long double error) {
        Point2D p(this->size, data);
        if (this->first == nullptr) {
            this->first = new Point2D(p.x, p.y);
        }
        else if (this->last == nullptr) {
            this->last = new Point2D(p.x, p.y);
        }

        if (this->last != nullptr) {
            Line n_line = this->__approximate(this->first, &p);
            if (this->__verify(this->l_cvx, this->u_cvx, n_line, error) > error) {
                return this->line->subs(this->size-1);
            }
            
            this->last->x = p.x; this->last->y = p.y;
            if (this->line != nullptr) delete this->line;
            this->line = new Line(n_line.get_slope(), n_line.get_intercept());    
        }

        this->size++;
        this->u_cvx.append(Point2D(p.x, p.y + error));
        this->l_cvx.append(Point2D(p.x, p.y - error));

        return INFINITY;
    }

    void InterpolateSegment::merge(Swab* neighbor) {
        InterpolateSegment* next = (InterpolateSegment*) neighbor;
        this->u_cvx.concat(next->u_cvx, this->size - 1);
        this->l_cvx.concat(next->l_cvx, this->size - 1);

        this->size += next->size - 1;
        this->last->x = this->size - 1;
        this->last->y = next->last->y;

        Line n_line = this->__approximate(this->first, this->last);
        if (this->line != nullptr) delete this->line;
        this->line = new Line(n_line.get_slope(), n_line.get_intercept());
    }

    long double InterpolateSegment::merge_cost(Swab* neighbor, long double error) {
        InterpolateSegment* next = (InterpolateSegment*) neighbor;
        // merge cost is related to the maximum individual error
        Line n_line = this->__approximate(this->first, next->last, this->size-1);
        long double error_1 = this->__verify(this->l_cvx, this->u_cvx, n_line, error);
        long double error_2 = this->__verify(next->l_cvx, next->u_cvx, n_line, error, this->size-1);

        return error_1 > error_2 ? error_1 : error_2;
    }
    // End: interpolate segment implementation

    // Begin: regression segment implementation
    Line RegressionSegment::__approximate(std::vector<long double>& data) {
        Eigen::MatrixXd A(data.size(), 2);
        Eigen::VectorXd b(data.size());

        for (int i = 0; i < data.size(); ++i) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
            b(i) = data[i];
            A(i, 0) = 1; A(i, 1) = i;
        }
        
        Eigen::VectorXd c = (((A.transpose() * A).inverse()) * A.transpose()) * b;
        return Line(c(1), c(0));
    }

    RegressionSegment::~RegressionSegment() {
        this->u_cvx.clear();
        this->l_cvx.clear();
        this->window.clear();
        if (this->line != nullptr) delete this->line;
    }

    BinObj* RegressionSegment::serialize(bool& first) {
        BinObj* obj = new BinObj;
        obj->put(VariableByteEncoding::encode(this->size));
        obj->put((float) this->line->get_intercept());
        obj->put((float) this->line->get_slope());

        return obj;
    }

    long double RegressionSegment::approximate(long double data, long double error) {
        Point2D p(this->size, data);
        this->window.push_back(p.y);

        if (this->window.size() == 2) {
            Line n_line = Line::line(
                Point2D(0, this->window[0]), Point2D(1, this->window[1])
            );

            this->line = new Line(n_line.get_slope(), n_line.get_intercept());    
        }
        else if (this->window.size() > 2) {
            Line n_line = this->__approximate(this->window);
            long double flag = this->__verify(this->l_cvx, this->u_cvx, n_line, error);
            if (
                flag > error || 
                std::abs(n_line.subs(p.x)-p.y) > error
            ) return this->line->subs(this->size-1);
            
            if (this->line != nullptr) delete this->line;
            this->line = new Line(n_line.get_slope(), n_line.get_intercept());    
        }

        this->size++;
        this->u_cvx.append(Point2D(p.x, p.y + error));
        this->l_cvx.append(Point2D(p.x, p.y - error));

        return INFINITY;
    }

    void RegressionSegment::merge(Swab* neighbor) {
        RegressionSegment* next = (RegressionSegment*) neighbor;
        this->window.insert(this->window.end(), next->window.begin(), next->window.end());
        this->u_cvx.concat(next->u_cvx, this->size);
        this->l_cvx.concat(next->l_cvx, this->size);

        this->size += next->size;
        Line n_line = this->__approximate(this->window);
        if (this->line != nullptr) delete this->line;
        this->line = new Line(n_line.get_slope(), n_line.get_intercept());
    }

    long double RegressionSegment::merge_cost(Swab* neighbor, long double error) {
        RegressionSegment* next = (RegressionSegment*) neighbor;
        std::vector<long double> s = this->window;
        s.insert(s.end(), next->window.begin(), next->window.end());
        // merge cost is related to the maximum individual error
        Line n_line = this->__approximate(s);

        long double error_1 = this->__verify(this->l_cvx, this->u_cvx, n_line, error);
        long double error_2 = this->__verify(next->l_cvx, next->u_cvx, n_line, error, this->size);

        return error_1 > error_2 ? error_1 : error_2;
    }
    // End: regression segment implementation

    // Begin: compression
    void Compression::__bottom_up() {
        std::vector<long double> m_err;
        for (int i=0; i<this->n_segment-1; i++) {
            Swab* current = this->segments[i];
            Swab* next = this->segments[i+1];
            m_err.push_back(current->merge_cost(next, this->error));
        }

        auto it = std::min_element(m_err.begin(), m_err.end());
        while (*it <= this->error && this->segments.size() > 1) {
            int index = std::distance(m_err.begin(), it);
            Swab* current = this->segments[index];
            Swab* next = this->segments[index+1];
            current->merge(next); delete next;
            
            m_err.erase(m_err.begin() + index);
            this->segments.erase(this->segments.begin() + index + 1);

            if (index != 0) m_err[index-1] = this->segments[index-1]->merge_cost(current, this->error);
            if (index < m_err.size()) m_err[index] = current->merge_cost(this->segments[index+1], this->error);
            it = std::min_element(m_err.begin(), m_err.end());  
        }
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
        this->n_segment = atoi(params[2]);

        if (this->mode == "interpolate") this->segment = new InterpolateSegment();
        else if (this->mode == "regression") this->segment = new RegressionSegment(); 
    }

    void Compression::finalize() {
        if (this->segment->line != nullptr) {
            this->segments.push_back(this->segment);
        }

        if (this->segments.size() > 1) {
            this->__bottom_up();
        }

        while (this->segments.size() != 0) {
            this->yield();
        }
        this->segments.clear();
    }

    void Compression::compress(Univariate* data) {
        long double tail = this->segment->approximate(data->get_value(), this->error);
        if (tail != INFINITY) {
            this->segments.push_back(this->segment);
            if (this->mode == "interpolate") {
                this->segment = new InterpolateSegment();
                this->segment->approximate(tail, this->error);
                this->segment->approximate(data->get_value(), this->error);
            }
            else if (this->mode == "regression") {
                this->segment = new RegressionSegment();
                this->segment->approximate(data->get_value(), this->error);
            } 
        }

        // Perform bottom up once n_segment sliding window segments have accumulated
        if (this->segments.size() == this->n_segment) {
            this->__bottom_up();
            this->yield();
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = this->segments[0]->serialize(this->first);

        delete this->segments[0];
        this->segments.erase(this->segments.begin());

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

        if (this->mode == "interpolate") {
            int start = 1;
            if (this->pivot == INFINITY) {
                start = 0;
                this->pivot = compress_data->getFloat();
            }

            unsigned long length = VariableByteEncoding::decode(compress_data);
            float last = compress_data->getFloat();
            Line line = Line::line(Point2D(0, this->pivot), Point2D(length-1, last));

            // std::cout << "decompress: " << length << " " << line.get_slope() << " " << line.get_intercept() << "\n";

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

            this->basetime += (length - 1) * this->interval;
            this->pivot = last;
        }
        else if (this->mode == "regression") {
            unsigned long length = VariableByteEncoding::decode(compress_data);
            float intercept = compress_data->getFloat();
            float slope = compress_data->getFloat();

            for (int i=0; i<length; i++) {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                    base_obj->pushData(std::to_string(slope * i + intercept));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime + i * this->interval));
                    obj->pushData(std::to_string(slope * i + intercept));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
            }

            this->basetime += length * this->interval;
        }
        
        return base_obj;
    }
    // End: decompression
};