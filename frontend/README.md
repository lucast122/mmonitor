# MMonitor Dashboard Frontend

Modern React + TypeScript dashboard for MMonitor metagenomics analysis.

## Tech Stack

- **React 18** - UI framework
- **TypeScript** - Type safety
- **Vite** - Build tool and dev server
- **TanStack Query** - Server state management
- **Zustand** - Client state management
- **Tailwind CSS** - Styling
- **Shadcn/ui** - UI components
- **Plotly.js** - Charting
- **React Router** - Routing

## Getting Started

### Prerequisites

- Node.js 18+
- npm or yarn

### Installation

```bash
# Install dependencies
npm install

# Start development server
npm run dev
```

The development server will start at http://localhost:5173

### Environment Variables

Copy `.env.example` to `.env` and configure:

```bash
cp .env.example .env
```

### Building for Production

```bash
npm run build
```

Built files will be in the `dist` directory.

## Project Structure

```
src/
├── api/           # API client and endpoints
├── components/
│   ├── ui/        # Base UI components (Shadcn)
│   ├── layout/    # Layout components (Sidebar, Header)
│   ├── filters/   # Filter components
│   └── charts/    # Chart components (Plotly)
├── hooks/         # Custom React hooks
├── lib/           # Utility functions
├── pages/         # Page components
├── stores/        # Zustand stores
└── types/         # TypeScript types
```

## Available Scripts

- `npm run dev` - Start development server
- `npm run build` - Build for production
- `npm run preview` - Preview production build
- `npm run lint` - Run ESLint

## API Integration

The frontend connects to the Django REST API at `/api/v1/`. In development, Vite proxies requests to `http://localhost:8000`.

## Features

- **Dashboard Overview** - Quick stats and navigation
- **Taxonomy Analysis** - Stacked bar, grouped bar, heatmap, area charts
- **Quality Control** - QC metrics tables and visualizations
- **Diversity Analysis** - Alpha and beta diversity metrics
- **MAG Analysis** - MAG quality scatter plots and details
