use std::ops::Index;

pub struct BitVector {
    data: Vec<u64>,
}

impl BitVector {
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub fn iter(&self) -> impl Iterator<Item = &u64> {
        self.data.iter()
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn count(&self) -> usize {
        self.data
            .iter()
            .fold(0, |acc, n| acc + n.count_ones() as usize)
    }
}

impl Default for BitVector {
    fn default() -> Self {
        Self::new()
    }
}

impl From<&[bool]> for BitVector {
    fn from(value: &[bool]) -> Self {
        let mut bv = BitVector::new();
        for (i, v) in value.iter().enumerate() {
            if bv.data.len() < i / 64 + 1 {
                bv.data.push(0);
            }
            bv.data[i / 64] |= (*v as u64) << (i % 64);
        }
        bv
    }
}

impl Index<usize> for BitVector {
    type Output = u64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl std::fmt::Debug for BitVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in &self.data {
            writeln!(f, "{:064b}", row)?;
        }
        Ok(())
    }
}
