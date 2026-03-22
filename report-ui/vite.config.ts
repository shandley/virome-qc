import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";
import { resolve } from "path";

export default defineConfig({
  plugins: [react(), tailwindcss()],
  resolve: {
    alias: {
      "@": resolve(__dirname, "./src"),
    },
  },
  build: {
    // Single file output for embedding in virome-qc
    rollupOptions: {
      output: {
        // Inline everything into a single JS file
        manualChunks: undefined,
        inlineDynamicImports: true,
        entryFileNames: "report.js",
        assetFileNames: "report.[ext]",
      },
    },
    // Inline assets under 100KB
    assetsInlineLimit: 100000,
    cssCodeSplit: false,
  },
});
