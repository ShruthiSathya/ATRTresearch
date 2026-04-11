import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  server: {
    port: 3000,
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
        // Strip the /api prefix before forwarding to FastAPI
        // FastAPI routes are /analyze, /generate_ai_analysis, etc. (no /api prefix)
        rewrite: (path) => path.replace(/^\/api/, '')
      }
    }
  }
})