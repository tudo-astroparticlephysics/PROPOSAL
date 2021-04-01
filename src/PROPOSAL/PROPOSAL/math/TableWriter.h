
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class Table_write {
   public:
    Table_write(std::string path, bool binary) : path_(path), binary_(binary) {
        if (binary_)
            file_.open(path_.append(".bin"), std::ios::out | std::ios::binary);
        else
            file_.open(path_.append(".txt"), std::ios::out);
    }

    template <typename T>
    void write(T data) {
        if (binary_)
            file_.write(reinterpret_cast<char*>(&data), sizeof(data));
        else
            file_ << data << '\n';
    }

    void close() { file_.close(); }

   private:
    std::string path_;
    bool binary_;
    std::fstream file_;
};

class Table_read {
   public:
    Table_read(std::string path, bool binary) : path_(path), binary_(binary) {
        if (binary_)
            file_.open(path_.append(".bin"), std::ios::in | std::ios::binary);
        else
            file_.open(path_.append(".txt"), std::ios::in);
        file_.seekg(file_.beg);
    }

    template <typename T>
    void read(T& data) {
        if (binary_)
            file_.read(reinterpret_cast<char*>(&data), sizeof(data));
        else {
            file_ >> data;
        }
    }

    void jump(int n) {
        for (int i = 0; i < n; ++i) {
            file_.ignore(256, '\n');
        }
    }

    void close() { file_.close(); }

   private:
    std::string path_;
    bool binary_;
    std::fstream file_;
    // int length_;
};
