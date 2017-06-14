#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <iterator>

#include "args.hxx"


namespace args
{

template<typename Out>
void split(const std::string &s, char delim, Out result)
{
  std::stringstream ss;
  ss.str(s);
  std::string item;

  while (std::getline(ss, item, delim)) {
    if (!item.empty()) {
      *(result++) = item;
    }
  }
}

std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

struct DoubleVectorReader
{
  void operator()(const std::string &name, const std::string &value, std::vector<double> &destination)
  {
    destination.clear();

    std::vector<std::string> elems = split(value, ' ');

    try {
      for (const std::string &elem : elems) {
        destination.push_back(std::stod(elem));
      }
    }
    catch (const std::invalid_argument &err) {
      throw args::ParseError(err.what());
    }
  }
};

}  // args
