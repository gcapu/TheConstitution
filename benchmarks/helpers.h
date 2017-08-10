#pragma once

#include <stdlib.h>
#include <limits.h>
#include <assert.h>


namespace TC
{

//Converts a string to int, float or double.
template <typename T> T convert(const std::string& val);

template <> int convert<int>(const std::string& val) {
  long temp = strtol(val.c_str(), NULL, 10);
  assert(temp <= INT_MAX);
  return static_cast<int>(temp) ;
  }
template <> float convert<float>(const std::string& val) {
  return strtof(val.c_str(), NULL);
  }
template <> double convert<double>(const std::string& val) {
  return strtod(val.c_str(), NULL);
  }

bool strStartsWith(const std::string& str, const std::string& start){
  if (start.size() <= str.size() && std::equal(start.begin(), start.end(), str.begin()))
    return true;
  else 
    return false;
  }

//Appends the tokens generated from line into the tokens vector.
//An optional number of values can be ignored.
//Returns true if the line ends with a separator, false otherwise.
template <typename T>
bool tokenize(std::vector<T>& tokens, const std::string& line, int ignore = 0, char separator = ','){
  std::stringstream lineStream(line);
  std::string cell;
  //ingnore the first few
  for(int i = 0; i<ignore; i++)
    std::getline(lineStream, cell, separator);
  //split the line into tokens
  while (std::getline(lineStream, cell, separator))
    tokens.push_back(convert<T>(cell));
  return !lineStream && cell.empty();
  }

} //namespace TC