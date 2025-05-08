// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>

#ifdef __cplusplus
  #include <yaml-cpp/yaml.h>

extern "C"
{
  typedef YAML::Node Yaml;
  typedef YAML::iterator YamlIterator;
#endif

  /// @brief Interoperatble string type
  struct StringT
  {
    char* ptr_;
    int size_;
  };

  /// @brief Interoperable array type for strings
  struct StringArrayT
  {
    StringT* ptr_;
    int size_;
  };

  /// @brief Interoperable array type for doubles
  struct DoubleArrayT
  {
    double* ptr_;
    int size_;
  };

  /// @brief Interoperable array type for YAML nodes
  struct NodeArrayT
  {
    Yaml** ptr_;
    int size_;
  };

  /// @brief Creates a YAML node from a string
  /// @param yaml_string YAML in string form
  /// @return pointer to the new YAML node
  Yaml* YamlCreateFromString(const char* yaml_string);

  /// @brief Creates a YAML node from a YAML file
  /// @param file_path path to the YAML file
  /// @return pointer to the new YAML node
  Yaml* YamlCreateFromFile(const char* file_path);

  /// @brief Outputs a YAML node to a file
  /// @param node YAML node to output
  /// @param file_path path to file to create (any existing file will be overwritten)
  void YamlToFile(Yaml* node, const char* file_path);

  /// @brief Returns the number of child elements in the node
  ///        This works for vectors and maps
  /// @param node YAML node to return size of
  /// @return number of node elements
  int YamlSize(Yaml* node);

  /// @brief Returns an iterator to the first child node
  /// @param node YAML node to iterate over
  /// @return beginning iterator
  YamlIterator* YamlBegin(Yaml* node);

  /// @brief Returns an iterator to one element past the last child node
  /// @param node YAML node to iterator over
  /// @return ending iterator
  YamlIterator* YamlEnd(Yaml* node);

  /// @brief Increments a YAML iterator
  /// @param iter YAML iterator to increment
  /// @param end YAML iterator one element past end
  /// @return true if incremented iter < end, false otherwise
  bool YamlIncrement(YamlIterator* iter, YamlIterator* end);

  /// @brief Checks if a YAML iterator is at the end
  /// @param iter YAML iterator to check
  /// @param end YAML iterator one element past end
  /// @return true if iter == end, false otherwise
  bool YamlAtEnd(YamlIterator* iter, YamlIterator* end);

  /// @brief Returns the key associated with a YAML iterator
  /// @param iter YAML iterator to return key for
  /// @return key as a c string
  StringT YamlKey(YamlIterator* iter);

  /// @brief Returns a sub-node
  /// @param node parent YAML node
  /// @param key key to find
  /// @param found true if successful, false otherwise
  /// @return sub-node
  Yaml* YamlGetNode(Yaml* node, const char* key, bool& found);

  /// @brief Gets a string from a YAML node
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return Pointer to string as const char array
  StringT YamlGetString(Yaml* node, const char* key, bool& found);

  /// @brief Gets an integer from a YAML node
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return integer value
  int YamlGetInt(Yaml* node, const char* key, bool& found);

  /// @brief Gets a float from a YAML node
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return float value
  float YamlGetFloat(Yaml* node, const char* key, bool& found);

  /// @brief Gets a double from a YAML node
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return double value
  double YamlGetDouble(Yaml* node, const char* key, bool& found);

  /// @brief Gets a boolean from a YAML node
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return boolean value
  bool YamlGetBool(Yaml* node, const char* key, bool& found);

  /// @brief Gets an array of strings from a YAML node
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return string array
  StringArrayT YamlGetStringArray(Yaml* node, const char* key, bool& found);

  /// @brief Gets an array of doubles from a YAML node
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return double array
  DoubleArrayT YamlGetDoubleArray(Yaml* node, const char* key, bool& found);

  /// @brief Gets an array of YAML nodes from a YAML node
  /// @details It is expected that the caller takes ownership of the individual
  ///          pointers to YAML nodes in the array
  /// @param node YAML node
  /// @param key key to search for
  /// @param found true if successful, false otherwise
  /// @return node array
  NodeArrayT YamlGetNodeArray(Yaml* node, const char* key, bool& found);

  /// @brief Gets a node from a YAML iterator
  /// @param iter YAML iterator
  /// @return YAML node
  Yaml* YamlGetNodeFromIterator(YamlIterator* iter);

  /// @brief Gets a string from a YAML iterator
  /// @param iter YAML iterator
  /// @return string as a c string
  StringT YamlGetStringFromIterator(YamlIterator* iter);

  /// @brief Gets an int from a YAML iterator
  /// @param iter YAML iterator
  /// @return integer value
  int YamlGetIntFromIterator(YamlIterator* iter);

  /// @brief Gets a float from a YAML iterator
  /// @param iter YAML iterator
  /// @return float value
  float YamlGetFloatFromIterator(YamlIterator* iter);

  /// @brief Gets a double from a YAML iterator
  /// @param iter YAML iterator
  /// @return double value
  double YamlGetDoubleFromIterator(YamlIterator* iter);

  /// @brief Gets a boolean from a YAML iterator
  /// @param iter YAML iterator
  /// @return boolean value
  bool YamlGetBoolFromIterator(YamlIterator* iter);

  /// @brief Gets an array of strings from a YAML iterator
  /// @param iter YAML iterator
  /// @return string array
  StringArrayT YamlGetStringArrayFromIterator(YamlIterator* iter);

  /// @brief Adds a YAML node to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value YAML node to add
  void YamlAddNode(Yaml* node, const char* key, Yaml* value);

  /// @brief Adds a string to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value string to add
  void YamlAddString(Yaml* node, const char* key, const char* value);

  /// @brief Adds an int to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value integer to add
  void YamlAddInt(Yaml* node, const char* key, int value);

  /// @brief Adds a float to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value float to add
  void YamlAddFloat(Yaml* node, const char* key, float value);

  /// @brief Adds a double to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value double to add
  void YamlAddDouble(Yaml* node, const char* key, double value);

  /// @brief Adds a boolean to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value boolean to add
  void YamlAddBool(Yaml* node, const char* key, bool value);

  /// @brief Adds an array of strings to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value string array to add
  void YamlAddStringArray(Yaml* node, const char* key, StringArrayT value);

  /// @brief Adds an array of doubles to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value double array to add
  void YamlAddDoubleArray(Yaml* node, const char* key, DoubleArrayT value);

  /// @brief Adds an array of YAML nodes to a YAML node
  /// @param node YAML node
  /// @param key key to apply value to
  /// @param value node array to add
  void YamlAddNodeArray(Yaml* node, const char* key, NodeArrayT value);

  /// @brief Copies a YAML node
  /// @param node YAML node to copy
  /// @return pointer to the new YAML node
  Yaml* YamlCopyNode(Yaml* node);

  /// @brief Copies a YAML node to a string
  /// @param node YAML node to copy
  /// @return pointer to the new string
  StringT YamlToString(Yaml* node);

  /// @brief Merges one YAML node into another
  /// @param dest destination YAML node
  /// @param src source YAML node
  /// @return true if successful, false otherwise
  bool YamlMergeNode(Yaml* dest, const Yaml* src);

  /// @brief Cleans up memory for a YAML node
  /// @param ptr Node pointer to free memory for
  void YamlDeleteNode(Yaml* ptr);

  /// @brief Cleans up memory for a char array
  /// @param string String to free memory for
  void YamlDeleteString(StringT string);

  /// @brief Cleans up memory for an array of strings
  /// @param array array to free memory for
  void YamlDeleteStringArray(StringArrayT array);

  /// @brief Cleans up memory for an array of doubles
  /// @param array array to free memory for
  void YamlDeleteDoubleArray(DoubleArrayT array);

  /// @brief Cleans up memory for an array of YAML nodes
  /// @details It is expected that the caller retains ownership of the
  ///          individual node pointers in the array
  /// @param array array to free memory for
  void YamlDeleteNodeArray(NodeArrayT array);

  /// @brief Cleans up memory for a YAML iterator
  /// @param ptr Iterator to free memory for
  void YamlDeleteIterator(YamlIterator* ptr);

#ifdef __cplusplus
}
#endif