pub struct PointMutation {
    pub chromosome: String,
    pub position: usize,
    pub reference: char,
    pub alternative: char,
}

pub struct Indel {
    pub chromosome: String,
    pub position: usize,
    pub anchor: char,
    pub deleted: Vec<char>,
    pub inserted: Vec<char>,
}

impl Indel {
    pub fn is_inframe(&self) -> bool {
        let net_length: isize = self.inserted.len() as isize - self.deleted.len() as isize;
        net_length.rem_euclid(3) == 0
    }
}
