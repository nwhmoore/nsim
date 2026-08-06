#[derive(Default)]
pub struct VectorSeries {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

#[derive(Debug, Default)]
pub struct KahanAccumulator {
    sum: f64,
    correction: f64,
}

impl KahanAccumulator {
    pub fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let next_sum = self.sum + adjusted;
        self.correction = (next_sum - self.sum) - adjusted;
        self.sum = next_sum;
    }

    #[must_use]
    pub fn total(&self) -> f64 {
        self.sum
    }
}
