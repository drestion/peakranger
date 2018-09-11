/*
 * profilezoom.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: xfeng
 */

#include "profilezoom.h"

namespace {
void get_sum(std::vector<double>& v,
             double& s) {
    for (uint32_t i = 0; i < v.size(); i++) {
        s += v[i];
    }
}
}

const
void profile_zoom::zoom_out(const uint32_t ws,
                            const data_type v,
                            data_type& r) {

//    rt_assert(ws)
//
//    /*
//     * get the (pos, r_count) pair
//     * is pos inside the previous window?
//     *          yes-> add to w_k r_count
//     *          no -> add to w_k+1 r_coount;
//     *                output w_k
//     */
//    uint32_t _w_b = 0;
//    uint32_t _w_e = _w_b + ws - 1;
//    uint32_t _w = 0;
//    uint32_t i = 0;
//    uint32_t _w_c = 0;
//    std::pair<uint32_t, int64_t> _e;
//    uint32_t _p = 0;
//    uint32_t _s = 0;
////    std::cout << "_w_b " << _w_b << " : _w_e " << _w_e << std::endl;
//    for (; i < v.size(); i++) {
//        //element_type& _e = v[i];
//        _p = v[i].first;
//        _s = v[i].getS();
//        _w = _p / ws;
////        std::cout << "_w_b " << _w_b << " : _w_e " << _w_e << std::endl;
//        if (_p >= _w_b && _p <= _w_e) {
//            _w_c += (uint32_t) _s;
////            std::cout << "_p " << _p << " : _w " << _w << std::endl;
////            std::cout << "adding" << std::endl;
//        } else {
////            std::cout << "_p " << _p << " : _w " << _w << std::endl;
////            std::cout << "pushing" << std::endl;
//
//            _e.first = (uint32_t) ((_w_b + _w_e) / 2);
//            _e.getS() = (int64_t) (_w_c);
//            r.push_back(_e);
//
//            _w_b = _w * ws;
//            _w_e = _w_b + ws - 1;
//            _w_c = _s;
////            std::cout << "adding" << std::endl;
//        }
//    }
//    _e.first = (uint32_t) ((_w_b + _w_e) / 2);
//    _e.getS() = (int64_t) (_w_c);
//    r.push_back(_e);
}

const void profile_zoom::zoom_out(const uint32_t ws,
                                  const uint32_t os,
                                  const data_type v,
                                  data_type & r) {

//    LOG_DEBUG1("Entering profile_zoom::zoom_out");
//
//    if (os >= ws || ws == 0 || v.size() < 1) {
//        return;
//    }
//
//    LOG_DEBUG2("starting\n");
//
//    uint32_t _w_b = 0;
//    uint32_t _w_e = _w_b + ws - 1;
//    uint32_t i = 0;
//    uint32_t j = 0;
//    element_type _e;
//    uint32_t _p = 0;
//    double _s = 0;
//    double _sab = 0;
//    double _san = 0;
//    double _sbm = 0;
//    double _snm = 0;
//    bool _last = false;
//    std::vector<uint32_t> _c;
//    std::queue<double> _vs;
//    std::queue<uint32_t> _vp;
//    std::vector<double> _van;
//    std::vector<double> _vbm;
//
//    //find the first points of all new clusters.
//    _c.push_back(v[0].first);
//    LOG_DEBUG2("cluster : "<<v[0].first<<"\n");
//    i = 1;
//    for (; i < v.size(); i++) {
//        _p = v[i].first;
//        if (_p > ws) {
//            if (_p - ws > _c.back() && _p - ws > v[i - 1].first) {
//                _c.push_back(_p);
//                LOG_DEBUG2( "cluster : "<<_p<<"\n");
//            }
//        }
//    }
//    ////find the next new cluster
//    //find the first point of a new cluster
//    //printout the first window
//    //slide to next window
//    //get the points inside this window
//    //printout this window
//    //slide to next window
//    //check if this window contains any points
//    //no points? exit
//
//    ////no more clusters? exit
//    i = 0;
//    j = 0;
//    for (; i < _c.size(); i++) {
//        //establish the firstwindow of this cluster
//        LOG_DEBUG2("got a new cluster:"<<_c[i]<<"\n");
//        _p = _c[i];
//        _w_e = _p;
//        if (_w_e + 1 < ws) {
//            _w_b = 0;
//        } else {
//            //the begining points
//            _w_b = _w_e - ws + 1;
//        }
//        if (_w_b == 0 && _w_e == 0) {
//            _w_b = 0;
//            _w_e = 1;
//        }LOG_DEBUG2("established first window "<<_w_b<<":"<<_w_e<<"\n");
//        while (true) {
//
//            //collect points in this window
//            for (; j < v.size(); j++) {
//
//                _p = v[j].first;
//                _s = v[j].getS();
//                LOG_DEBUG2("  got "<<j<<"th point "<<_p<<":"<<_s<<"->");
//                if (_p >= _w_b && _p <= _w_e) {
//                    _vs.push(_s);
//                    _vp.push(_p);
//                    _vbm.push_back(_s);
//                    LOG_DEBUG2("collected new points in this window : "<<_p <<":"<<_s<<"\n");
//                } else {
//                    LOG_DEBUG2("this point passed the current window\n ");
//                    break;
//                }
//            }
//            if (_vs.size() > 0) {
//
//                //window not empty
//                get_sum(_vbm,
//                        _sbm);
//                get_sum(_van,
//                        _san);
//                _snm = _sab - _san + _sbm;
//                LOG_DEBUG2(" _sab:"<<_sab<<", _san:"<<_san<<", _sbm:"<<_sbm<<" _snm = _sab -_san+_sbm _ "<<_snm<<"\n");
//                _e.first = _w_e;
//                _e.getS() = _snm / ws;
//                LOG_DEBUG2(" this window is not empty. Got _e:"<<_e.first <<":"<<_e.getS()<<"\n");
//                //print out this window
//
//                r.push_back(_e);
//                if (_last) {
//                    break;
//                }
//                { //move to next window
//                    _van.clear();
//                    _vbm.clear();
//                    _sab = _snm;
//                    _san = 0;
//                    _sbm = 0;
//                    _w_b = _w_b + (ws - os);
//                    _w_e = _w_b + ws - 1;
//                    LOG_DEBUG2("moved to a new window "<<_w_b<<":"<<_w_e<<"\n");
//                    if (_w_e >= v.back().first) {
//                        LOG_DEBUG2("********this window's boundary reaches or passes the last data point. \n ");
//                        _last = true;
//                    }rt_assert(_vs.size() == _vp.size())
//                    while (_vp.size()) {
//                        LOG_DEBUG2("Checking old point in the previous window: "<<_vp.front()<<"->");
//                        if (_vp.front() < _w_b) {
//                            LOG_DEBUG2("popped out.\n");
//                            _van.push_back(_vs.front());
//                            _vs.pop();
//                            _vp.pop();
//
//                        } else {
//                            LOG_DEBUG2("in it. no more popout.\n");
//                            break;
//                        }
//                    }
//                }
//            } else {
//                //window empty, move to next cluster
//                LOG_DEBUG2("moved to next cluster.\n");
//                break;
//            }
//        }
//    }rt_assert(j == v.size()) LOG_DEBUG1("QUIT profile_zoom::zoom_out");
}

const void profile_zoom::smooth(const uint32_t ws,
                                const uint32_t os,
                                const data_type v,
                                data_type & r) {

    LOG_DEBUG1("Entering profile_zoom::zoom_out");

    if (os >= ws || ws == 0 || v.size() < 1) {
        LOG_DEBUG1("INPUT INVALID profile_zoom::zoom_out");
        return;
    }
    LOG_DEBUG2("window size "<<ws<<"\n");
    LOG_DEBUG2("overlap size "<<os<<"\n");
    LOG_DEBUG2("total wig points : "<<v.size()<<"\n");
    LOG_DEBUG2("the fisrt wig point : "<<v.front().getP()<<"\n");
    LOG_DEBUG2("the last wig point : "<<v.back().getP()<<"\n");
    uint32_t _w_b = 0;
    uint32_t _hw = ws / 2;
    uint32_t _w_e = _w_b + ws - 1;
    uint32_t i = 0;
    uint32_t j = 0;
    element_type _e;
    uint32_t _p = 0;
    double _s = 0;
    double _sab = 0;
    double _san = 0;
    double _sbm = 0;
    double _snm = 0;
    bool _last = false;
    std::vector<uint32_t> _c;
    std::queue<double> _vs;
    std::queue<uint32_t> _vp;
    std::vector<double> _van;
    std::vector<double> _vbm;

    //find the first points of all new clusters.
    _c.push_back(v[0].getP());
    LOG_DEBUG2("cluster : "<<v[0].getP()<<"\n");
    i = 1;
    for (; i < v.size(); i++) {
        _p = v[i].getP();
        if (_p > ws) {
            if (_p - ((ws / 2) + 1) > _c.back()
            && _p - ((ws / 2) + 1) > v[i - 1].getP()) {
                _c.push_back(_p);
                LOG_DEBUG2( "cluster : "<<_p<<"\n");
            }
        }
    }
    ////find the next new cluster
    //find the first point of a new cluster
    //printout the first window
    //slide to next window
    //get the points inside this window
    //printout this window
    //slide to next window
    //check if this window contains any points
    //no points? exit

    ////no more clusters? exit
    i = 0;
    j = 0;
    for (; i < _c.size(); i++) {
        //establish the first window of this cluster
        LOG_DEBUG2("got a new cluster:"<<_c[i]<<"\n");
        _p = _c[i];
        _w_e = _p + _hw;
        if (_w_e + 1 < ws) {
            _w_b = 0;
        } else {
            //the beginning points
            _w_b = _w_e - ws + 1;
        }
        if (_w_b == 0 && _w_e == 0) {
            _w_b = 0;
            _w_e = 1;
        }

        LOG_DEBUG2("established first window "<<_w_b<<":"<<_w_e<<"\n");

        while (true) {

            //collect points in this window
            for (; j < v.size(); j++) {

                _p = v[j].getP();
                _s = v[j].getS();
                LOG_DEBUG2("  got the "<<j+1<<"th point "<<_p<<":"<<_s<<"->");
                if (_p >= _w_b && _p <= _w_e) {
                    _vs.push(_s);
                    _vp.push(_p);
                    _vbm.push_back(_s);
                    LOG_DEBUG2("collected new points in this window : "<<_p <<":"<<_s<<"\n");
                } else {
                    LOG_DEBUG2("this point passed the current window\n ");
#ifdef USE_LOGGING
                    if(_vs.size() < 1) {
                        LOG_DEBUG2("*****This window passed "
                        "without collecting any points\n");
                    }
#endif
                    break;
                }
            }
            if (_vs.size() > 0) {

                //window not empty
                get_sum(_vbm,
                        _sbm);
                get_sum(_van,
                        _san);
                _snm = _sab - _san + _sbm;
                LOG_DEBUG2(" _sab:"<<_sab<<", _san:"<<_san<<", _sbm:"<<_sbm<<" _snm = _sab -_san+_sbm _ "<<_snm<<"\n");
                _e.setP((_w_e + _w_b) / 2);
                _e.setS(_snm / ws);
                LOG_DEBUG2(" this window is not empty. Got _e:"<<_e.getP() <<":"<<_e.getS()<<"\n");
                //print out this window

                r.push_back(_e);
//                if (_last) {
//                    LOG_DEBUG2(" Finished processing the last window.\n");
//                    break;
//                }
                { //move to next window
                    _van.clear();
                    _vbm.clear();
                    _sab = _snm;
                    _san = 0;
                    _sbm = 0;
                    _w_b = _w_b + (ws - os);
                    _w_e = _w_b + ws - 1;
                    LOG_DEBUG2("moved to a new window "<<_w_b<<":"<<_w_e<<"\n");
                    if (_w_b > v.back().getP()) {

                        LOG_DEBUG2("********this window's boundary reaches"
                        " or passes the last data point. \n ");

                        _last = true;
                    }

                    assert(_vs.size() == _vp.size());

                    while (_vp.size()) {
                        LOG_DEBUG2("Checking old point in the previous window");
                        if (_vp.front() < _w_b) {
                            LOG_DEBUG2("The front point "<<_vp.front()
                            <<" is not in the window. "
                            "popped it out.\n");
                            _van.push_back(_vs.front());
                            _vs.pop();
                            _vp.pop();

                        } else {
                            LOG_DEBUG2("The front point is in the window. no more popout.\n");
                            break;
                        }
                    }
                }
            } else {
#ifdef USE_LOGGING
                //window empty, move to next cluster
                LOG_DEBUG2("moved to next cluster.\n");
                if(i == (_c.size()-1)) {
                    LOG_DEBUG2("hooray!!! this is the last cluster!");
                }
#endif
                break;
            }
        }
    }

    //todo: w=10k u =1k triggers this
    rt_assert(j == v.size())

    LOG_DEBUG1("QUIT profile_zoom::smooth");
}
