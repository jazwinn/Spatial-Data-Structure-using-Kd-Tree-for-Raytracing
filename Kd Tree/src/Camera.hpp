#ifndef CAMERA_HPP
#define CAMERA_HPP

#include "Math.hpp"

namespace CS350 {
    class Camera
    {
      private:
        glm::mat4  mView{};
        glm::mat4  mProj{};
        glm::mat4  mProjView{};
        glm::vec3  mPosition{};
        glm::vec3  mTarget{};
        glm::vec3  mUp{};
        float      mFovDeg{};
        glm::ivec2 mSize{};
        float      mNear{};
        float      mFar{};
        bool       mIsOrthogonal = false;

      public:
        Camera();
        void Update();
        void SetProjection(float fov_deg, glm::ivec2 window_size, float near, float far);
        void SetProjection(glm::ivec2 window_size, float near, float far);
        void SetPosition(glm::vec3 position);
        void SetTarget(glm::vec3 position);
        void SetUp(glm::vec3 up);

        glm::vec3        position() const { return mPosition; }
        glm::vec3        target() const { return mTarget; }
        glm::vec3        up() const { return mUp; }
        glm::mat4 const& view() const { return mView; }
        glm::mat4&       view() { return mView; }
        glm::mat4 const& proj() const { return mProj; }
        glm::mat4&       proj() { return mProj; }
        glm::mat4 const& viewProj() const { return mProjView; }
        void             set_size(glm::ivec2 v) { mSize = v; }
        glm::ivec2       size() const { return mSize; }
        void             set_near(float v) { mNear = v; }
        float            near() const { return mNear; }
        void             set_far(float v) { mFar = v; }
        float            far() const { return mFar; }
        float            fov_deg() const { return mFovDeg; }
        void             set_fov_deg(float f) { mFovDeg = f; }

      private:
        void EnsureOrthonormal();
    };
}

#endif // CAMERA_HPP
