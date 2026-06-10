#include <cmath>
#include <iostream>

template <typename T>
class Vector3
{
public:
    T       x;
    T       y;
    T       z;
    double coords[3];
    Vector3() = default;

    Vector3(T x, T y, T z) : x(x), y(y), z(z), coords{ x,y,z }{}


    double Magnitude() {
        return sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
    }

    Vector3<T> Cross(const Vector3<T>& other) {
        return Vector3<T>(this->y * other.z - this->z * other.y, this->z * other.x - this->x * other.z, this->x * other.y - this->y * other.x);
    }

    Vector3<T> Copy() {
        return Vector3<T>(this->x, this->y, this->z);
    }


    void Set(T x, T y, T z) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->coords[0] = x;
        this->coords[1] = y;
        this->coords[2] = z;
    }

    friend Vector3<T> operator*(const Vector3<T>& lhs, T r) {
        return Vector3<T>(lhs.x * r, lhs.y * r, lhs.z * r);
    }
    friend Vector3<T> operator*(T r, const Vector3<T>& rhs) {
        return Vector3<T>(rhs.x * r, rhs.y * r, rhs.z * r);
    }

    friend Vector3<T> operator+(const Vector3<T>& lhs, const Vector3<T>& rhs) {
        return Vector3<T>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    friend Vector3<T> operator-(const Vector3<T>& lhs, const Vector3<T>& rhs) {
        return Vector3<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }
    friend Vector3<T> operator/(const Vector3<T>& lhs, const Vector3<T>& rhs) {
        return Vector3<T>(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z);
    }
    friend Vector3<T> operator/(const Vector3<T>& lhs, double den) {
        return Vector3<T>(lhs.x / den, lhs.y / den, lhs.z / den);
    }
    friend std::ostream& operator<<(std::ostream& os, const Vector3& v) {
        os << "Vector3: " << v.x << ", " << v.y << ", " << v.z;
        return os;
    }
};

struct Camera
{
    double              near, far;
    double              angle;
    Vector3<double>     position;
    Vector3<double>     focus;
    Vector3<double>     up;
};

Camera SetupCamera(void);

struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char* colors;  // size is 3*numBins
    double* opacities; // size is numBins

    // Take in a value and applies the transfer function.
    // Step #1: figure out which bin "value" lies in.
    // If "min" is 2 and "max" is 4, and there are 10 bins, then
    //   bin 0 = 2->2.2
    //   bin 1 = 2.2->2.4
    //   bin 2 = 2.4->2.6
    //   bin 3 = 2.6->2.8
    //   bin 4 = 2.8->3.0
    //   bin 5 = 3.0->3.2
    //   bin 6 = 3.2->3.4
    //   bin 7 = 3.4->3.6
    //   bin 8 = 3.6->3.8
    //   bin 9 = 3.8->4.0
    // and, for example, a "value" of 3.15 would return the color in bin 5
    // and the opacity at "opacities[5]".

    int GetBin(double value) {
        double delta = (this->max - this->min) / (double) this->numBins;

        for (int i = 0; i < this->numBins; i++) {
            double cv = this->min + delta * i;
            double nv = this->min + delta * (i + 1);

            if (value >= cv && value < nv) {
                return i;
            }
        }
        return -1;
    }

    static int GetBin(double value, double max, double min, double n_bins) {
        double delta = (max - min) / (double)n_bins;

        for (int i = 0; i < n_bins; i++) {
            double cv = min + delta * i;
            double nv = min + delta * (i + 1);

            if (value >= cv && value < nv) {
                return i;
            }
        }
        return -1;
    }

    void ApplyTransferFunction(double value, unsigned char* RGB, double& opacity)
    {
        int bin = GetBin(value);
        if (bin == -1) {
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
            opacity = 0;
            return;
        }
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
    }

    static void ApplyTransferFunction(double value, unsigned char* RGB, unsigned char* colors, int bin)
    {
        if (bin == -1) {
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
            return;
        }
        RGB[0] = colors[3 * bin + 0];
        RGB[1] = colors[3 * bin + 1];
        RGB[2] = colors[3 * bin + 2];
        //opacity = opacities[bin];
    }
};



TransferFunction SetupTransferFunction(void);
