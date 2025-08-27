#include "Logging.hpp"
#include "Shapes.hpp"

namespace {
    template <typename T>
    std::istream& GenericVecRead(std::istream& is, T& v) {
        for (int i = 0; i < T::length(); ++i) {
            is >> v[i];
            if (i + 1 != T::length()) {
                is.ignore(1, ',');
            }
        }
        return is;
    }

    template <typename T>
    std::ostream& GenericVecWrite(std::ostream& os, T const& v) {
        for (int i = 0; i < T::length(); ++i) {
            os << v[i];
            if (i + 1 != T::length()) {
                os << ", ";
            }
        }
        return os;
    }
}

namespace glm {
    std::istream& operator>>(std::istream& is, vec2& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, vec2 const& v) { return GenericVecWrite(os, v); }
    std::istream& operator>>(std::istream& is, vec3& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, vec3 const& v) { return GenericVecWrite(os, v); }
    std::istream& operator>>(std::istream& is, vec4& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, vec4 const& v) { return GenericVecWrite(os, v); }
    std::istream& operator>>(std::istream& is, ivec2& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, ivec2 const& v) { return GenericVecWrite(os, v); }
    std::istream& operator>>(std::istream& is, ivec3& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, ivec3 const& v) { return GenericVecWrite(os, v); }
    std::istream& operator>>(std::istream& is, ivec4& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, ivec4 const& v) { return GenericVecWrite(os, v); }
    std::istream& operator>>(std::istream& is, mat3& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, mat3 const& v) { return GenericVecWrite(os, v); }
    std::istream& operator>>(std::istream& is, mat4& v) { return GenericVecRead(is, v); }
    std::ostream& operator<<(std::ostream& os, mat4 const& v) { return GenericVecWrite(os, v); }
}

namespace CS350 {
    std::ostream& operator<<(std::ostream& os, Ray const& ray) {
        os << fmt::format("ray.start={{{}}}; ray.dir={{{}}};", ray.start, ray.dir);
        return os;
    }
    std::ostream& operator<<(std::ostream& os, Triangle const& tri) {
        os << fmt::format("A:({}), B({}), C:({})", tri[0], tri[1], tri[2]);
        return os;
    }
    std::ostream& operator<<(std::ostream& os, Aabb const& aabb) {
        os << fmt::format("min:({}), max({})", aabb.min, aabb.max);
        return os;
    }
}
