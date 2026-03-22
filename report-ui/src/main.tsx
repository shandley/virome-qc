import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import "./index.css";
import { App } from "./App";

// Load passport data from embedded script tag
const dataEl = document.getElementById("passport-data");
const passportData = dataEl ? JSON.parse(dataEl.textContent || "{}") : null;

// Detect if this is a batch report (array) or single sample
const isBatch = Array.isArray(passportData);

// Apply saved theme
const saved = localStorage.getItem("virome-qc-theme");
const prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
if (saved === "dark" || (!saved && prefersDark)) {
  document.documentElement.classList.add("dark");
}

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <App data={passportData} isBatch={isBatch} />
  </StrictMode>
);
