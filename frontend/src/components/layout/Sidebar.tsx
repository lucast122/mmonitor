import { NavLink } from 'react-router-dom'
import {
  LayoutDashboard,
  Dna,
  BarChart3,
  Microscope,
  HardDrive,
  Layers,
  GitCompare,
} from 'lucide-react'
import { cn } from '@/lib/utils'

const navItems = [
  { to: '/', icon: LayoutDashboard, label: 'Overview' },
  { to: '/taxonomy', icon: Dna, label: 'Taxonomy' },
  { to: '/qc', icon: BarChart3, label: 'Quality Control' },
  { to: '/diversity', icon: Layers, label: 'Diversity' },
  { to: '/correlations', icon: GitCompare, label: 'Correlations' },
  { to: '/mags', icon: Microscope, label: 'MAGs' },
]

export function Sidebar() {
  return (
    <aside className="hidden w-64 flex-shrink-0 border-r bg-card lg:block">
      <div className="flex h-16 items-center border-b px-6">
        <HardDrive className="mr-2 h-6 w-6 text-primary" />
        <span className="text-xl font-bold">MMonitor</span>
      </div>

      <nav className="space-y-1 p-4">
        {navItems.map((item) => (
          <NavLink
            key={item.to}
            to={item.to}
            end={item.to === '/'}
            className={({ isActive }) =>
              cn(
                'flex items-center gap-3 rounded-lg px-3 py-2 text-sm font-medium transition-colors',
                isActive
                  ? 'bg-primary text-primary-foreground'
                  : 'text-muted-foreground hover:bg-muted hover:text-foreground'
              )
            }
          >
            <item.icon className="h-5 w-5" />
            {item.label}
          </NavLink>
        ))}
      </nav>
    </aside>
  )
}
