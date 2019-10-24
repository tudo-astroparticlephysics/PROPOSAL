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
