#include "piecewise-approximation/linear.hpp"

namespace IOnlyPLA {
    // Begin: compression
    Eigen::VectorXd Compression::__feasible_rectangle() {
        Eigen::VectorXd x(3); Eigen::VectorXd c(3);                                       
        Eigen::VectorXd b(this->u_cvx.size()+this->l_cvx.size());  
        Eigen::MatrixXd A(this->u_cvx.size()+this->l_cvx.size(), 3);                    

        c(0) = 0; c(1) = 1; c(2) = -1;
        for (int i=0; i<this->u_cvx.size(); i++) {
            A(i, 0) = -this->u_cvx.at(i).x; 
            A(i, 1) = -1; 
            A(i, 2) = 0; 
            b(i) = -this->u_cvx.at(i).y;
        }

        for (int i=0; i<this->l_cvx.size(); i++) {
            A(i+this->u_cvx.size(), 0) = this->l_cvx.at(i).x; 
            A(i+this->u_cvx.size(), 1) = 0; 
            A(i+this->u_cvx.size(), 2) = 1; 
            b(i+this->u_cvx.size()) = this->l_cvx.at(i).y;
        }
        
        double minobj = this->sdlp.linprog(c, A, b, x);
        if (minobj == INFINITY || minobj == -INFINITY) {
            std::cout << "wtf\n";
        }

        return x;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
    }

    void Compression::finalize() {
        if (this->length >= 2) this->yield();
        
        if (this->pivot != nullptr) delete this->pivot;
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        this->u_cvx.clear();
        this->l_cvx.clear();
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        Eigen::VectorXd c = this->__feasible_rectangle();
        Line l_bound(c(0), c(1));
        Line u_bound(c(0), c(2));

        long u_right = static_cast<long>(std::floor(u_bound.subs(this->length)));
        long l_right = static_cast<long>(std::ceil(l_bound.subs(this->length)));
        long u_left = static_cast<long>(std::floor(u_bound.get_intercept()));
        long l_left = static_cast<long>(std::ceil(l_bound.get_intercept()));

        if (u_left >= l_left && u_right >= l_right) {
            long embedded = this->length << 1 | 0;
            obj->put(VariableByteEncoding::encode(embedded));
            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode(l_left)));
            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode(u_right)));
        }
        else {
            long scale = std::ceil(1/(c(2) - c(1)));
            u_right = static_cast<long>(std::floor(scale*u_bound.subs(this->length)));
            l_left = static_cast<long>(std::ceil(scale*l_bound.get_intercept()));
            long embedded = this->length << 1 | 1;

            obj->put(VariableByteEncoding::encode(embedded));
            obj->put(VariableByteEncoding::encode(scale));
            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode(l_left)));
            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode(u_right)));
        }

        return obj;
    }
    
    void Compression::compress(Univariate* data) {
        Point2D p(this->length++, data->get_value());

        if (this->pivot == nullptr) {
            this->pivot = new Point2D(p.x, p.y);
            this->u_cvx.append(Point2D(p.x, p.y - this->error));
            this->l_cvx.append(Point2D(p.x, p.y + this->error));
        }
        else if (this->u_line == nullptr) {
            Line u_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y-this->error), 
                Point2D(p.x, p.y+this->error)
            );
            Line l_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y+this->error), 
                Point2D(p.x, p.y-this->error)
            );

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
            this->u_cvx.append(Point2D(p.x, p.y - this->error));
            this->l_cvx.append(Point2D(p.x, p.y + this->error));
        }
        else {
            bool complete = this->l_line->subs(p.x) - p.y - this->error > 0.0000001 
                || p.y - this->error - this->u_line->subs(p.x) > 0.0000001;
            
            if (!complete) {
                bool update_u = p.y + this->error < this->u_line->subs(p.x);
                bool update_l = p.y - this->error > this->l_line->subs(p.x);

                int u_index = -1;
                int l_index = -1;
                Line* prev_u_line = this->u_line;
                Line* prev_l_line = this->l_line;
                
                if (update_u) {
                    long double min_slp = INFINITY;

                    for (int i=0; i<this->u_cvx.size(); i++) {
                        Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y + this->error));
                        if (line.get_slope() < min_slp) {
                            min_slp = line.get_slope();
                            u_index = i;
                        }
                    }
                    Line line = Line::line(this->u_cvx.at(u_index), Point2D(p.x, p.y + this->error));
                    this->u_line = new Line(line.get_slope(), line.get_intercept());
                }
                if (update_l) {
                    long double max_slp = -INFINITY;

                    for (int i=0; i<this->l_cvx.size(); i++) {
                        Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y - this->error));
                        if (line.get_slope() > max_slp) {
                            max_slp = line.get_slope();
                            l_index = i;
                        }
                    }
                    Line line = Line::line(this->l_cvx.at(l_index), Point2D(p.x, p.y - this->error));
                    this->l_line = new Line(line.get_slope(), line.get_intercept());
                }
                
                bool identical = std::abs(this->u_line->get_slope() - this->l_line->get_slope()) < 0.0000001
                    && std::abs(this->u_line->get_intercept() - this->l_line->get_intercept()) < 0.0000001;

                if (identical) {
                    complete = true;
                    if (update_u) { delete this->u_line; this->u_line = prev_u_line; }
                    if (update_l) { delete this->l_line; this->l_line = prev_l_line; }
                }
                else {
                    if (update_u) { 
                        delete prev_u_line; 
                        this->u_cvx.erase_from_begin(u_index);
                        this->l_cvx.append(Point2D(p.x, p.y + this->error)); 
                    }
                    if (update_l) { 
                        delete prev_l_line; 
                        this->l_cvx.erase_from_begin(l_index);
                        this->u_cvx.append(Point2D(p.x, p.y - this->error)); 
                    }
                }
            }

            if (complete) {
                this->length--;
                this->yield();
                
                delete this->u_line; this->u_line = nullptr;
                delete this->l_line; this->l_line = nullptr;
                this->u_cvx.clear(); this->l_cvx.clear();

                this->pivot->x = 0; this->pivot->y = p.y;
                u_cvx.append(Point2D(0, p.y - this->error));
                l_cvx.append(Point2D(0, p.y + this->error));
                this->length = 1;
            }
        }
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        // Do nothing
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        unsigned long embedded = VariableByteEncoding::decode(compress_data);
        int flag = embedded & 1;
        int length = embedded >> 1;
        
        int scale = 1;
        if (flag) scale = VariableByteEncoding::decode(compress_data);
        long value_1 = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));
        long value_2 = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));

        Line line = Line::line(Point2D(0, (long double) value_1 / scale), 
            Point2D(length, (long double) value_2 / scale));

        for (int i=0; i<length; i++) {
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

        this->basetime += length * this->interval;
        return base_obj;
    }
    // End: decompression
};