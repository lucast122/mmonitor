import { type ClassValue, clsx } from 'clsx'
import { twMerge } from 'tailwind-merge'

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs))
}

export function formatNumber(num: number, decimals = 2): string {
  if (num >= 1e9) return (num / 1e9).toFixed(decimals) + 'B'
  if (num >= 1e6) return (num / 1e6).toFixed(decimals) + 'M'
  if (num >= 1e3) return (num / 1e3).toFixed(decimals) + 'K'
  return num.toFixed(decimals)
}

export function formatPercentage(num: number, decimals = 1): string {
  return num.toFixed(decimals) + '%'
}

export function downloadAsCSV(data: Record<string, unknown>[], filename: string) {
  if (data.length === 0) return

  const headers = Object.keys(data[0])
  const csvContent = [
    headers.join(','),
    ...data.map((row) =>
      headers
        .map((header) => {
          const value = row[header]
          // Escape quotes and wrap in quotes if contains comma
          const strValue = String(value ?? '')
          return strValue.includes(',') || strValue.includes('"')
            ? `"${strValue.replace(/"/g, '""')}"`
            : strValue
        })
        .join(',')
    ),
  ].join('\n')

  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' })
  const link = document.createElement('a')
  link.href = URL.createObjectURL(blob)
  link.download = filename
  link.click()
}

export function getColorPalette(count: number): string[] {
  // Extended color palette for charts
  const palette = [
    '#2563eb', // blue-600
    '#dc2626', // red-600
    '#16a34a', // green-600
    '#ca8a04', // yellow-600
    '#9333ea', // purple-600
    '#ea580c', // orange-600
    '#0891b2', // cyan-600
    '#db2777', // pink-600
    '#4f46e5', // indigo-600
    '#65a30d', // lime-600
    '#0d9488', // teal-600
    '#7c3aed', // violet-600
    '#f97316', // orange-500
    '#06b6d4', // cyan-500
    '#ec4899', // pink-500
    '#8b5cf6', // violet-500
    '#84cc16', // lime-500
    '#14b8a6', // teal-500
    '#f59e0b', // amber-500
    '#6366f1', // indigo-500
  ]

  if (count <= palette.length) {
    return palette.slice(0, count)
  }

  // Generate additional colors if needed
  const additionalColors: string[] = []
  for (let i = palette.length; i < count; i++) {
    const hue = (i * 137.5) % 360
    additionalColors.push(`hsl(${hue}, 70%, 50%)`)
  }

  return [...palette, ...additionalColors]
}
