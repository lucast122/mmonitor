import axios, { AxiosError, InternalAxiosRequestConfig } from 'axios'
import { useAuthStore } from '@/stores/authStore'

const API_BASE_URL = import.meta.env.VITE_API_URL || '/api/v1'

export const apiClient = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
})

// Request interceptor to add auth token
apiClient.interceptors.request.use(
  (config: InternalAxiosRequestConfig) => {
    const token = useAuthStore.getState().accessToken
    if (token) {
      config.headers.Authorization = `Bearer ${token}`
    }
    return config
  },
  (error) => Promise.reject(error)
)

// Response interceptor to handle token refresh
apiClient.interceptors.response.use(
  (response) => response,
  async (error: AxiosError) => {
    const originalRequest = error.config as InternalAxiosRequestConfig & { _retry?: boolean }

    // If 401 and not already retrying, try to refresh token
    if (error.response?.status === 401 && !originalRequest._retry) {
      originalRequest._retry = true

      const refreshToken = useAuthStore.getState().refreshToken

      if (refreshToken) {
        try {
          const response = await axios.post(`${API_BASE_URL}/auth/refresh/`, {
            refresh: refreshToken,
          })

          const { access } = response.data
          useAuthStore.getState().setAccessToken(access)

          // Retry original request with new token
          originalRequest.headers.Authorization = `Bearer ${access}`
          return apiClient(originalRequest)
        } catch (refreshError) {
          // Refresh failed, logout user
          useAuthStore.getState().logout()
          window.location.href = '/dashboard/login'
          return Promise.reject(refreshError)
        }
      } else {
        // No refresh token, logout user
        useAuthStore.getState().logout()
        window.location.href = '/dashboard/login'
      }
    }

    return Promise.reject(error)
  }
)

export default apiClient
