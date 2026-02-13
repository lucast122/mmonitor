import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  base: '/dashboard/',
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
      // Polyfill Node.js built-ins used by plotly.js transitive dependencies
      buffer: 'buffer/',
      stream: 'stream-browserify',
      assert: 'assert/',
    },
  },
  server: {
    port: 5173,
    host: '0.0.0.0',  // Listen on all interfaces
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
      },
    },
  },
  build: {
    outDir: 'dist',
    sourcemap: true,
  },
})
