// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/util/config_yaml.h>

#include <yaml-cpp/yaml.h>

#include <cstring>
#include <fstream>

Yaml* YamlCreateFromString(const char* yaml_string)
{
  return new YAML::Node(YAML::Load(yaml_string));
}

Yaml* YamlCreateFromFile(const char* file_path)
{
  return new YAML::Node(YAML::LoadFile(file_path));
}

void YamlToFile(Yaml* node, const char* file_path)
{
  std::ofstream file(file_path, std::ofstream::trunc);
  file << *node;
  file.close();
}

int YamlSize(Yaml* node)
{
  return node->size();
}

YamlIterator* YamlBegin(Yaml* node)
{
  return new YAML::iterator(node->begin());
}

YamlIterator* YamlEnd(Yaml* node)
{
  return new YAML::iterator(node->end());
}

bool YamlAtEnd(YamlIterator* iter, YamlIterator* end)
{
  return *iter == *end;
}

bool YamlIncrement(YamlIterator* iter, YamlIterator* end)
{
  return ++(*iter) != *end;
}

StringT YamlKey(YamlIterator* iter)
{
  StringT string;
  std::string str = (*iter)->first.as<std::string>();
  string.size_ = str.length();
  string.ptr_ = new char[string.size_ + 1];
  strcpy(string.ptr_, str.c_str());
  return string;
}

Yaml* YamlGetNode(Yaml* node, const char* key, bool& found)
{
  YAML::Node subnode = (*node)[key];
  found = subnode.IsDefined() && !subnode.IsScalar();
  return new YAML::Node(subnode);
}

StringT YamlGetString(Yaml* node, const char* key, bool& found)
{
  found = (*node)[key].IsDefined();
  StringT string;
  if (found)
  {
    std::string str = (*node)[key].as<std::string>();
    string.size_ = str.length();
    string.ptr_ = new char[string.size_ + 1];
    strcpy(string.ptr_, str.c_str());
    return string;
  }
  string.ptr_ = nullptr;
  string.size_ = 0;
  return string;
}

int YamlGetInt(Yaml* node, const char* key, bool& found)
{
  found = (*node)[key].IsDefined();
  if (found)
    return (*node)[key].as<int>();
  return 0;
}

float YamlGetFloat(Yaml* node, const char* key, bool& found)
{
  found = (*node)[key].IsDefined();
  if (found)
    return (*node)[key].as<float>();
  return 0.0f;
}

double YamlGetDouble(Yaml* node, const char* key, bool& found)
{
  found = (*node)[key].IsDefined();
  if (found)
    return (*node)[key].as<double>();
  return 0.0;
}

bool YamlGetBool(Yaml* node, const char* key, bool& found)
{
  found = (*node)[key].IsDefined();
  if (found)
    return (*node)[key].as<bool>();
  return false;
}

StringArrayT YamlGetStringArray(Yaml* node, const char* key, bool& found)
{
  StringArrayT array;
  array.size_ = 0;
  array.ptr_ = nullptr;
  YAML::Node array_node = (*node)[key];
  found = array_node.IsDefined();
  if (!found)
    return array;
  array.size_ = array_node.size();
  array.ptr_ = new StringT[array.size_];
  for (std::size_t i = 0; i < array_node.size(); ++i)
  {
    std::string str = array_node[i].as<std::string>();
    array.ptr_[i].size_ = str.length();
    array.ptr_[i].ptr_ = new char[str.length() + 1];
    strcpy(array.ptr_[i].ptr_, str.c_str());
  }
  return array;
}

DoubleArrayT YamlGetDoubleArray(Yaml* node, const char* key, bool& found)
{
  DoubleArrayT array;
  array.size_ = 0;
  array.ptr_ = nullptr;
  YAML::Node array_node = (*node)[key];
  found = array_node.IsDefined();
  if (!found)
    return array;
  array.size_ = array_node.size();
  array.ptr_ = new double[array.size_];
  for (std::size_t i = 0; i < array_node.size(); ++i)
  {
    array.ptr_[i] = array_node[i].as<double>();
  }
  return array;
}

NodeArrayT YamlGetNodeArray(Yaml* node, const char* key, bool& found)
{
  NodeArrayT array;
  array.size_ = 0;
  array.ptr_ = nullptr;
  YAML::Node array_node = (*node)[key];
  found = array_node.IsDefined();
  if (!found)
    return array;
  array.size_ = array_node.size();
  array.ptr_ = new YAML::Node*[array.size_];
  for (std::size_t i = 0; i < array_node.size(); ++i)
  {
    array.ptr_[i] = new YAML::Node(array_node[i].as<YAML::Node>());
  }
  return array;
}

Yaml* YamlGetNodeFromIterator(YamlIterator* iter)
{
  return (*iter)->IsDefined() ? new YAML::Node((*iter)->as<YAML::Node>()) : new YAML::Node((*iter)->second.as<YAML::Node>());
}

StringT YamlGetStringFromIterator(YamlIterator* iter)
{
  StringT string;
  std::string str = (*iter)->IsDefined() ? (*iter)->as<std::string>() : (*iter)->second.as<std::string>();
  string.size_ = str.length();
  string.ptr_ = new char[string.size_ + 1];
  strcpy(string.ptr_, str.c_str());
  return string;
}

int YamlGetIntFromIterator(YamlIterator* iter)
{
  return (*iter)->IsDefined() ? (*iter)->as<int>() : (*iter)->second.as<int>();
}

float YamlGetFloatFromIterator(YamlIterator* iter)
{
  return (*iter)->IsDefined() ? (*iter)->as<float>() : (*iter)->second.as<float>();
}

double YamlGetDoubleFromIterator(YamlIterator* iter)
{
  return (*iter)->IsDefined() ? (*iter)->as<double>() : (*iter)->second.as<double>();
}

bool YamlGetBoolFromIterator(YamlIterator* iter)
{
  return (*iter)->IsDefined() ? (*iter)->as<bool>() : (*iter)->second.as<bool>();
}

StringArrayT YamlGetStringArrayFromIterator(YamlIterator* iter)
{
  StringArrayT array;
  YAML::Node array_node = (*iter)->IsDefined() ? (*iter)->as<YAML::Node>() : (*iter)->second.as<YAML::Node>();
  array.size_ = array_node.size();
  array.ptr_ = new StringT[array.size_];
  for (std::size_t i = 0; i < array_node.size(); ++i)
  {
    std::string str = array_node[i].as<std::string>();
    array.ptr_[i].size_ = str.length();
    array.ptr_[i].ptr_ = new char[str.length() + 1];
    strcpy(array.ptr_[i].ptr_, str.c_str());
  }
  return array;
}

void YamlAddNode(Yaml* node, const char* key, Yaml* value)
{
  (*node)[key] = YAML::Clone(*value);
}

void YamlAddString(Yaml* node, const char* key, const char* value)
{
  (*node)[key] = value;
}

void YamlAddInt(Yaml* node, const char* key, int value)
{
  (*node)[key] = value;
}

void YamlAddFloat(Yaml* node, const char* key, float value)
{
  (*node)[key] = value;
}

void YamlAddDouble(Yaml* node, const char* key, double value)
{
  (*node)[key] = value;
}

void YamlAddBool(Yaml* node, const char* key, bool value)
{
  (*node)[key] = value;
}

void YamlAddStringArray(Yaml* node, const char* key, StringArrayT value)
{
  YAML::Node array;
  for (std::size_t i = 0; i < value.size_; ++i)
  {
    array.push_back(value.ptr_[i].ptr_);
  }
  (*node)[key] = array;
}

void YamlAddDoubleArray(Yaml* node, const char* key, DoubleArrayT value)
{
  YAML::Node array;
  for (std::size_t i = 0; i < value.size_; ++i)
  {
    array.push_back(value.ptr_[i]);
  }
  (*node)[key] = array;
}

void YamlAddNodeArray(Yaml* node, const char* key, NodeArrayT value)
{
  YAML::Node array;
  for (std::size_t i = 0; i < value.size_; ++i)
  {
    array.push_back(*(value.ptr_[i]));
  }
  (*node)[key] = array;
}

Yaml* YamlCopyNode(Yaml* node)
{
  return new YAML::Node(YAML::Clone(*node));
}

StringT YamlToString(Yaml* node)
{
  StringT string;
  YAML::Emitter out;
  out << *node;
  string.size_ = out.size();
  string.ptr_ = new char[string.size_ + 1];
  strcpy(string.ptr_, out.c_str());
  return string;
}

bool YamlMergeNode(Yaml* node, const Yaml* other)
{
  if (!node->IsMap() || !other->IsMap())
    return false;
  for (YAML::const_iterator it = (*other).begin(); it != (*other).end(); ++it)
  {
    std::string key = it->first.as<std::string>();
    if ((*node)[key].IsDefined() && (*node)[key].IsMap() && it->second.IsMap())
    {
      Yaml subnode = (*node)[key];
      if (!YamlMergeNode(&subnode, &it->second))
        return false;
      (*node)[key] = subnode;
    }
    else
    {
      if ((*node)[key].IsDefined() && !(*node)[key].is(it->second))
      {
        return false;
      }
      (*node)[key] = it->second;
    }
  }
  return true;
}

void YamlDeleteNode(Yaml* ptr)
{
  delete ptr;
}

void YamlDeleteString(StringT string)
{
  delete[] string.ptr_;
}

void YamlDeleteStringArray(StringArrayT array)
{
  if (!array.ptr_)
    return;
  for (std::size_t i = 0; i < array.size_; ++i)
  {
    delete[] array.ptr_[i].ptr_;
  }
  delete[] array.ptr_;
}

void YamlDeleteDoubleArray(DoubleArrayT array)
{
  delete[] array.ptr_;
}

void YamlDeleteNodeArray(NodeArrayT array)
{
  delete[] array.ptr_;
}

void YamlDeleteIterator(YamlIterator* ptr)
{
  delete ptr;
}