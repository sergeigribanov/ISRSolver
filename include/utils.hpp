#ifndef __UTILS_HPP__
#define __UTILS_HPP__
#include <TFile.h>

#include <iostream>
#include <string>

template <class T>
T* find_object(TFile* fl, const std::string& object_name) {
  T* object = dynamic_cast<T*>(fl->Get(object_name.c_str()));
  if (!object) {
    std::cerr << "[!] Object with name \"" << object_name << "\" of class "
              << T::Class_Name() << " is not found in file " << fl->GetName()
              << std::endl;
    exit(1);
  }
  return object;
}

#endif
