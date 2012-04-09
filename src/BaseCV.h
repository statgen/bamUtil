/*
 *  Copyright (C) 2010-2012  Christian Fuchsberger,
 *                           Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __Base_CV_H__
#define __Base_CV_H__

#include <stdint.h>
#include <string>
#include <vector>


// will be refactored!

struct baseCV {
    uint8_t q;
    uint16_t pos;
    uint8_t tile;
    uint8_t read;
    uint8_t prebase;
    uint8_t nexbase;
    uint8_t curbase;
    uint8_t obs;
    uint8_t ref;
    std::string rg; //readgroup
    std::vector<uint16_t> covariates;
    std::vector<uint8_t> bitsize;

    void init(std::string readGroup, uint16_t k0,uint16_t k1,uint16_t k2, uint16_t k3, uint8_t k4, uint8_t k5, uint8_t k6,uint8_t k7, uint8_t k8)
    {
     rg = readGroup;
     q = k0;
     pos= k1;
     tile = k2;
     read= k3;
     prebase= k4;
     nexbase= k5;
     curbase= k6;
     obs = k7;
     ref = k8;
     covariates.push_back(k0);
     covariates.push_back(k1);
     //covariates.push_back(k2);
     covariates.push_back(k3);
     covariates.push_back(k4);
     //covariates.push_back(k5);
     covariates.push_back(k6);
    }
    void setCovariates(){
        covariates.clear();
    	covariates.push_back(q);
        covariates.push_back(pos);
        //covariates.push_back(tile);
        covariates.push_back(read);
        covariates.push_back(prebase);
        //covariates.push_back(nexbase);
        covariates.push_back(curbase);
    }
};

#endif
