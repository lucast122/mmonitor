declare module 'd3-horizon' {
  interface HorizonOptions {
    useCanvas?: boolean
  }

  export default class Horizon {
    constructor(element: HTMLElement, options?: HorizonOptions)
    width(px: number): this
    height(px: number): this
    data(data: [number, number][]): this
    x(accessor: ((d: unknown) => number) | string): this
    y(accessor: ((d: unknown) => number) | string): this
    xMin(value: number): this
    xMax(value: number): this
    yExtent(value: number): this
    yScaleExp(value: number): this
    yAggregation(fn: (vals: number[]) => number): this
    bands(count: number): this
    mode(mode: 'offset' | 'mirror'): this
    positiveColors(colors: string[]): this
    positiveColorStops(stops: number[]): this
    negativeColors(colors: string[]): this
    negativeColorStops(stops: number[]): this
    interpolationCurve(curve: unknown): this
    duration(ms: number): this
    tooltipContent(fn: (d: { x: number; y: number; points?: unknown[] }) => string | null): this
    onHover(fn: (d: { x: number; y: number; points?: unknown[] } | null) => void): this
    onClick(fn: (d: { x: number; y: number; points?: unknown[] } | null) => void): this
  }
}
