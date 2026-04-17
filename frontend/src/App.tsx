import { Routes, Route, Navigate } from 'react-router-dom'
import { useAuthStore } from '@/stores/authStore'
import { Layout } from '@/components/layout/Layout'
import { LoginPage } from '@/pages/LoginPage'
import { DashboardPage } from '@/pages/DashboardPage'
import { TaxonomyPage } from '@/pages/TaxonomyPage'
import { QCPage } from '@/pages/QCPage'
import { DiversityPage } from '@/pages/DiversityPage'
import { CorrelationsPage } from '@/pages/CorrelationsPage'
import { MAGsPage } from '@/pages/MAGsPage'

function ProtectedRoute({ children }: { children: React.ReactNode }) {
  const isAuthenticated = useAuthStore((state) => state.isAuthenticated)

  if (!isAuthenticated) {
    return <Navigate to="/login" replace />
  }

  return <>{children}</>
}

function App() {
  return (
    <Routes>
      <Route path="/login" element={<LoginPage />} />
      <Route
        path="/"
        element={
          <ProtectedRoute>
            <Layout />
          </ProtectedRoute>
        }
      >
        <Route index element={<DashboardPage />} />
        <Route path="taxonomy" element={<TaxonomyPage />} />
        <Route path="qc" element={<QCPage />} />
        <Route path="diversity" element={<DiversityPage />} />
        <Route path="correlations" element={<CorrelationsPage />} />
        <Route path="mags" element={<MAGsPage />} />
      </Route>
    </Routes>
  )
}

export default App
