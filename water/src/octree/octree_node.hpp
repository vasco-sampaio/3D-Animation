#pragma once

#include "cgp/cgp.hpp"

template <typename T>
class Octree;

template <typename Container>
concept ContainerConcept = requires(Container c) {
    typename Container::value_type;
    typename Container::iterator;
};

template <typename T>
class OctreeNode {
private:
    OctreeNode const* _parent = nullptr; // avoid using const_cast that removes const qualifier from an object
    mutable OctreeNode* _children[8] = { nullptr }; // property that can be modified even in a const object
    unsigned int _code;
    mutable T _value;

    const OctreeNode* searchValue(unsigned int mortonCode) const {
        if (mortonCode == _code) {
            return this;
        }

        int index = mortonCode & 0x7;

        if (!_children[index]) {
            _children[index] = new OctreeNode(mortonCode);
            _children[index]->_parent = this;
        }

        return _children[index]->searchValue(mortonCode);
    }

    OctreeNode(unsigned int code) : _code(code), _value(T{}) {}
    OctreeNode(unsigned int code, const T& value) : _code(code), _value(value) {}

    ~OctreeNode() {
        for (int i = 0; i < 8; ++i) {
            delete _children[i];
        }
    }

public:
    const OctreeNode* get_parent() const { return _parent; }
    const OctreeNode* get_children() const { return _children; }
    T& get_value() const { return _value; }

    std::enable_if_t <ContainerConcept<T>, T> get_siblings_values() const {
        T siblingsValues;

        for (int i = 0; i < 8; ++i) {
            if (this->_parent->_children[i] != this) {
                const auto& siblingValue = this->_parent->_children[i]->get_value();
                siblingsValues.insert(siblingValue.begin(), siblingValue.end());
            }
        }

        return siblingsValues;
    }

    friend class Octree<T>;
};
