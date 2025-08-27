#include "Camera.hpp"
#include "Math.hpp"
#include <glm/ext/matrix_clip_space.hpp>

namespace CS350 {
    Camera::Camera()
    {
        SetPosition({ 10, 10, 10 });
        SetTarget({});
        SetProjection(60.0f, { 4, 3 }, 0.1f, 1000.0f);
        Update();
    }

    void Camera::Update()
    {
        mView = glm::lookAt(mPosition, mTarget, mUp);

        // Orthogonal
        if (mIsOrthogonal) {
            mProj = glm::ortho(float(mSize.x) * -0.5f, float(mSize.x) * 0.5f, float(mSize.y) * -0.5f, float(mSize.y) * 0.5f, mNear, mFar);
        }
        // Perspective
        else {
            float aspect_ratio = static_cast<float>(mSize.x) / static_cast<float>(mSize.y);
            mProj              = glm::perspective(glm::radians(mFovDeg), aspect_ratio, mNear, mFar);
        }

        mProjView = mProj * mView;
    }

    void Camera::SetProjection(float fov_deg, glm::ivec2 window_size, float near, float far)
    {
        mFovDeg       = fov_deg >= 1 ? fov_deg : 60.0f;
        mSize         = window_size.x > 0 ? window_size : glm::ivec2{ 640, 480 };
        mNear         = near >= 0.01f ? near : 0.01f;
        mFar          = far > near ? far : near + 100.0f;
        mIsOrthogonal = false;
    }

    void Camera::SetProjection(glm::ivec2 window_size, float near, float far)
    {
        mIsOrthogonal = true;
        mSize         = window_size.x > 0 ? window_size : glm::ivec2{ 640, 480 };
        mNear         = near >= 0.01f ? near : 0.01f;
        mFar          = far > near ? far : near + 100.0f;
    }

    void Camera::SetPosition(glm::vec3 position)
    {
        mPosition = position;
        EnsureOrthonormal();
    }

    void Camera::SetTarget(glm::vec3 position)
    {
        mTarget = position;
        EnsureOrthonormal();
    }

    void Camera::SetUp(glm::vec3 up)
    {
        mUp = up;
        EnsureOrthonormal();
    }

    void Camera::EnsureOrthonormal()
    {
        glm::vec3 v    = mTarget - mPosition;
        glm::vec3 side = glm::cross(v, { 0, 1, 0 });
        if (glm::all(glm::epsilonEqual(side, glm::vec3{ 0 }, 0.01f))) {
            mUp = { 1, 0, 0 };
        } else {
            mUp = { 0, 1, 0 };
        }
    }
}
