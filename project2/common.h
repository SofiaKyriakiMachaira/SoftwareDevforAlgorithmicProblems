
#if !defined(MyHeader_H)
#define MyHeader_H
using namespace std;
// header should be defined just once
struct POINTDISTANCE { // Saves position of point and its euclidian distance
  int point;
  double distance;
};

struct PDCompareFunctionDistance { // Compares distance between 2 points
    inline bool operator() (const POINTDISTANCE& pd1, const POINTDISTANCE& pd2){
        return (pd1.distance < pd2.distance);
    }
};

struct PDCompareFunctionPoint { // Compares position between 2 points
    inline bool operator() (const POINTDISTANCE& pd1, const POINTDISTANCE& pd2){
        return (pd1.point < pd2.point);
    }
};

struct PDCompareFunctionEqualPoint { // Checks equality between two points of POINTDISTANCE
  explicit PDCompareFunctionEqualPoint(const int& s) : point(s) {}
  inline bool operator() (const POINTDISTANCE& pd1){
    return (pd1.point == point);
  }
  int point;
};

#endif