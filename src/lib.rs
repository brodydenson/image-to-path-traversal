use std::fmt;
use std::rc::Rc;
use std::collections::LinkedList;
use image::{ImageBuffer, Luma};
use std::path::Path;
type CellRef = Rc<Cell>;

pub struct Canvas {
    grid: Grid,
    skeleton: Grid,
    path: LinkedList<CellRef>,
    dx_dt_path: LinkedList<i8>,
    dy_dt_path: LinkedList<i8>,
}

impl Canvas {
    pub fn new() -> Canvas {
        Canvas { grid: Grid::new(), skeleton: Grid::new(), path: LinkedList::new(), dx_dt_path: LinkedList::new(), dy_dt_path: LinkedList::new() }
    }

    pub fn build_from_img(&mut self, str_src_path: &str) {
        println!("Reading image...");

        let src_path = Path::new(str_src_path);

        let img = image::open(src_path).expect("Error opening input image").to_luma8();

        let (width, height) = img.dimensions();
        self.grid.set_dim((width, height));
        println!("Images size (width, height): {}, {}", width, height);

        for (x, y, pixel) in img.enumerate_pixels() {
            if (x != 0 && y != 0)
                && (x != width-1 && y != height-1)
                && pixel.0[0] < 128 {
                self.grid.set_cell(x, y, true);
            }
        }

        // println!("Building outline...");
        // self.build_outline();
        println!("Building skeleton... (may take a while)");
        self.build_skeleton();
        println!("Building path... (may take a while)");
        self.build_path();
    }

    fn build_img(&self, grid: &Grid, path: &Path) -> Result<(), image::ImageError> {
        let (width, height) = self.grid.dim;
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            let in_grid = grid.at(x, y);
            *pixel = match in_grid {
                Some(_) => Luma([0]),       
                _ => Luma([225]),
            }
        }

        return new_img.save(path);
    }

    pub fn save(&self) {
        let original_dir = Path::new("build/original.png");
        let skeleton_dir = Path::new("build/skeleton.png");
        let path_dir = Path::new("build/path.png");

        // Saving original image
        println!("Saving original image...");
        self.build_img(&self.grid, original_dir).expect("Error saving original image");

        // Saving skeleton image
        println!("Saving skeleton image...");
        self.build_img(&self.skeleton, skeleton_dir).expect("Error saving skeleton image");

        // Saving path image
        println!("Saving path image...");
        let (width, height) = self.grid.dim;
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            let in_path = self.path.iter().position(|cell| cell.x == x && cell.y == y);

            *pixel = match in_path {
                Some(i) => Luma([(200.0 * (i as f64 / self.path.len() as f64)) as u8]),       
                _ => Luma([225]),
            }
        }
        new_img
            .save(path_dir)
            .expect("Error saving path image");
    }

    // fn build_outline(&mut self) {
    //     self.outline = self.grid.select(|cell| self.grid.get_adjacent_cross(&cell).len() < 4);
    // }

    fn build_skeleton(&mut self) {
        self.skeleton = self.grid.clone();

        let mut prev_remove_count = 0;
        let mut curr_remove_count = self.thin();
        while curr_remove_count != prev_remove_count {
            prev_remove_count = curr_remove_count;
            curr_remove_count = self.thin();
            println!("Remove Count: {}", curr_remove_count);
            println!("Skeleton Length: {}", self.skeleton.len);
            continue;
        }
    }

    // Zhang Suen thinning algorithm (https://arxiv.org/pdf/1710.03025.pdf)
    fn thin(&mut self) -> u32 {
        let test_1 = |grid: &Grid, cell: CellRef| -> bool {
            // (0) The pixel is black and has eight neighbours
            // (1) 2 <= B(P1) <= 6 
            // (2) A(P1) = 1
            // (3) At least one of P2 and P4 and P6 is white
            // (4) At least one of P4 and P6 and P8 is white
            let (x, y) = cell.get_pos();
            // (0)
            if !(x + 1 < grid.dim.0 && y + 1 < grid.dim.1 && x as i32 - 1 >= 0 && y as i32 - 1 >= 0) { return false }
            // (1)
            let taken_adjacent_n = grid.get_adjacent(&cell).len();
            if !(2 <= taken_adjacent_n && taken_adjacent_n <= 6) { return false }
            // (2)
            if !(grid.get_transitions(&cell) == 1) { return false }
            // (3)
            if !(grid.at(x, y + 1).is_none() ||
                grid.at(x + 1, y).is_none() ||
                grid.at(x, y - 1).is_none()) { return false }
            // (4)
            if !(grid.at(x + 1, y).is_none() ||
                grid.at(x, y - 1).is_none() ||
                grid.at(x - 1, y).is_none()) { return false }

            true
        };
        let test_2 = |grid: &Grid, cell: CellRef| -> bool {
            // (0) The pixel is black and has eight neighbours
            // (1) 2 <= B(P1) <= 6 
            // (2) A(P1) = 1
            // (3) At least one of P2 and P4 and P8 is white
            // (4) At least one of P2 and P6 and P8 is white
            let (x, y) = cell.get_pos();
            // (0)
            if !(x + 1 < grid.dim.0 && y + 1 < grid.dim.1 && x as i32 - 1 >= 0 && y as i32 - 1 >= 0) { return false }
            // (1)
            let taken_adjacent_n = grid.get_adjacent(&cell).len();
            if !(2 <= taken_adjacent_n && taken_adjacent_n <= 6) { return false }
            // (2)
            if !(grid.get_transitions(&cell) == 1) { return false }
            // (3)
            if !(grid.at(x, y + 1).is_none() ||
                grid.at(x + 1, y).is_none() ||
                grid.at(x - 1, y).is_none()) { return false }
            // (4)
            if !(grid.at(x, y + 1).is_none() ||
                grid.at(x, y - 1).is_none() ||
                grid.at(x - 1, y).is_none()) { return false }

            true
        };

        let mut remove_count = 0;

        let viewed_grid = self.skeleton.clone();
        let selected = self.skeleton.select(|cell: &CellRef| test_1(&viewed_grid, Rc::clone(&cell)));
        for cell in selected {
            self.skeleton.set_cell(cell.x, cell.y, false);
            remove_count += 1;
        }

        let viewed_grid = self.skeleton.clone();
        let selected = self.skeleton.select(|cell: &CellRef| test_2(&viewed_grid, Rc::clone(&cell)));
        for cell in selected {
            self.skeleton.set_cell(cell.x, cell.y, false);
            remove_count += 1;
        }

        remove_count
    }

    fn build_path(&mut self) {
        let mut path: LinkedList<CellRef> = LinkedList::new();
        let mut prev_path_len = 0;
        let mut starting_cells = self.skeleton.select(|cell| self.skeleton.get_adjacent(cell).len() <= 2 && !path.contains(cell));
        starting_cells.append(&mut self.skeleton.select(|cell| self.skeleton.get_adjacent(cell).len() == 3 && !path.contains(cell)));
        while (path.len() > prev_path_len || prev_path_len == 0) && starting_cells.len() > 0 {
            println!("Path Length: {}", path.len());
            prev_path_len = path.len();
            let mut max_path = self.max_path(Rc::clone(&starting_cells[0]), &LinkedList::new(), &path, 0);
            path.append(&mut max_path);
            starting_cells = self.skeleton.select(|cell| self.skeleton.get_adjacent(cell).len() <= 2 && !path.contains(cell));
            starting_cells.append(&mut self.skeleton.select(|cell| self.skeleton.get_adjacent(cell).len() == 3 && !path.contains(cell)));
        }
        self.path = path;

        let mut dx_dt_path: LinkedList<i8> = LinkedList::new();
        let mut dy_dt_path: LinkedList<i8> = LinkedList::new();

        fn distribute_add(list: &mut LinkedList<i8>, num: i32) {
            if num == 0 { list.push_back(0); return; }
            let val = if num > 0 { 1 } else { -1 };
            for _ in 0..num.abs() {
                list.push_back(val);
            }
        }
        // Basic constant derivative (mostly discontinuous)
        let mut prev_cell: &CellRef = self.path.front().unwrap();
        for cell in &self.path {
            if cell == prev_cell { continue; }
            distribute_add(&mut dx_dt_path, cell.x as i32 - prev_cell.x as i32);
            distribute_add(&mut dy_dt_path, cell.y as i32 - prev_cell.y as i32);
            prev_cell = cell;
        }

        self.dx_dt_path = dx_dt_path;
        self.dy_dt_path = dy_dt_path;
    }

    fn max_path(&self, cell: CellRef, taken: &LinkedList<CellRef>, past_taken: &LinkedList<CellRef>, depth: u8) -> LinkedList<CellRef> {
        if depth >= 8 { return taken.clone(); }

        let mut curr_taken = taken.clone();
        let mut adjacent: Vec<CellRef>;
        let mut curr_cell = Rc::clone(&cell);
        let get_taken_adjacent = |grid: &Grid, taken: &LinkedList<CellRef>, curr_cell: &CellRef| {
            grid.get_adjacent(curr_cell).iter().filter_map(|cell| {
                if cell != curr_cell && !taken.contains(cell) {
                    return Some(Rc::clone(cell));
                }
                None
            }).collect()
        };

        'travel: loop {
            adjacent = get_taken_adjacent(&self.skeleton, &curr_taken, &curr_cell);
            curr_taken.push_back(Rc::clone(&curr_cell));
            if adjacent.len() > 1 {
                // WTF IS THIS, I JUST CAME BACK TO THIS AND IM ALREADY CONSFUSED
                'outer: for cell_a in &adjacent {
                    let cell_a_adjacent = get_taken_adjacent(&self.skeleton, &curr_taken, cell_a);
                    for cell_b in &adjacent {
                        if cell_a != cell_b && !cell_a_adjacent.contains(&cell_b) { continue 'outer }
                    }
                    // Can't intersect more than one point
                    if past_taken.contains(&curr_cell) && past_taken.contains(&cell_a) { 
                        return curr_taken
                    }
                    curr_cell = Rc::clone(&cell_a);
                    continue 'travel
                }
                break 'travel
            }
            // Dead end
            if adjacent.len() == 0 { 
                return curr_taken
            }

            curr_cell = Rc::clone(&adjacent[0]);
        }


        // Splitting to multiple paths
        let mut max_taken: LinkedList<CellRef> = LinkedList::new();
        for cell in adjacent {
            // New path
            let new_taken = self.max_path(cell, &curr_taken, past_taken, depth+1);
            if new_taken.len() > max_taken.len() {
                max_taken = new_taken;
            }
        }
        
        max_taken
    }

    pub fn print(&self) {
        self.grid.print();
    }
}

#[derive(Eq, PartialEq, Clone)]
pub struct Grid {
    dim: (u32, u32),
    len: u32,
    cells: Vec<Vec<Option<CellRef>>>,
}
 
impl Grid {
    pub fn new() -> Grid {
        Grid { dim: (0, 0), cells: Vec::new(), len: 0 }
    }

    pub fn set_dim(&mut self, dim: (u32, u32)) {
        self.dim = dim;
        self.cells = vec![vec![None; dim.1 as usize]; dim.0 as usize];
    }

    pub fn set_cell(&mut self, row: u32, col: u32, some: bool) {
        let cell_opt = if some {
            if self.at(row, col).is_none() { self.len += 1 }
            Some(Rc::new(Cell { x: row as u32, y: col as u32 }))
        } else {
            if self.at(row, col).is_some() { self.len -= 1 }
            None
        };

        self.cells[row as usize][col as usize] = cell_opt;
    }

    pub fn select<F>(&self, f: F) -> Vec<CellRef> where
        F: Fn(&CellRef) -> bool {

        let mut selected_cells: Vec<CellRef> = Vec::new();
        let (height, width) = self.dim;

        for row in 0..height {
            for col in 0..width {
                if let Some(taken_cell) = self.at(row as u32, col as u32) {
                    if f(&taken_cell) {
                        selected_cells.push(taken_cell);
                    }
                }
            }
        }

        selected_cells
    }

    pub fn print(&self) {
        let (height, width) = self.dim;

        for row in 0..height {
            for col in 0..width {
                print!("{}",
                    match self.at(row as u32, col as u32) {
                        Some(_) => '#',
                        None => '.'
                    }
                )
            }
            print!("\n");
        }
    }

    pub fn print_vertices(&self, vertices: &Vec<CellRef>) {
        let (height, width) = self.dim;

        let find = |x: u32, y: u32| {
            for v in vertices {
                if v.get_pos() == (x, y) {
                    return true;
                }
            }
            false
        };

        for row in 0..height {
            for col in 0..width {
                if find(row as u32, col as u32) {
                    print!("#");
                } else {
                    print!(".");
                }
            }
            print!("\n");
        }
    }

    pub fn at(&self, row: u32, col: u32) -> Option<CellRef> {
        // Out of bounds
        if row >= self.dim.0 || col >= self.dim.1 { return None }

        let cell = &self.cells[row as usize][col as usize];
        match cell {
            Some(taken) => Some(Rc::clone(taken)),
            _ => None,
        }

    }

    pub fn get_adjacent(&self, cell: &Cell) -> Vec<CellRef> {
        let mut adjacent: Vec<CellRef> = Vec::new();
        let (x, y) = cell.get_pos();

        for row_offset in -1..=1 {
            for col_offset in -1..=1 {
                let row = x as i32 + row_offset;
                let col = y as i32 + col_offset;

                if row < 0 || col < 0 { continue }

                if let Some(adj_cell) = self.at(row as u32, col as u32) {
                    adjacent.push(adj_cell);
                }
            }
        }

        adjacent
    }

    pub fn get_transitions(&self, cell: &Cell) -> u8 {
        let mut none_to_sum_n = 0;
        let mut curr_some = true; // cell must be some
        let (x, y) = cell.get_pos();

        let circular_adjacent = vec![(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1), (0, 1)];
        for (x_offset, y_offset) in circular_adjacent {
            let cell = self.at((x as i32 + x_offset) as u32, (y as i32 + y_offset) as u32);
            if !curr_some && cell.is_some() {
                none_to_sum_n += 1;
            }
            curr_some = cell.is_some();
        }

        none_to_sum_n
    }

    pub fn get_adjacent_cross(&self, cell: &Cell) -> Vec<CellRef> {
        let mut adjacent: Vec<CellRef> = Vec::new();
        let (x, y) = cell.get_pos();

        for row_offset in -1..=1 {
            for col_offset in -1..=1 {
                // Adjecent not including the element itself, and diaginal
                if (row_offset != 0 && col_offset != 0)
                    || (row_offset == 0 && col_offset == 0) {
                    continue;
                }

                let row = x as i32 + row_offset;
                let col = y as i32 + col_offset;

                // Out of bounds
                if row < 0 || col < 0 {
                    continue;
                }

                let (row, col) = (row as usize, col as usize);

                // Out of bounds
                if row as u32 >= self.dim.0 || col as u32 >= self.dim.1 {
                    continue;
                }

                if let Some(taken_cell) = self.at(row as u32, col as u32) {
                    adjacent.push(taken_cell);
                }
            }
        }

        adjacent
    }
}

#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub struct Cell {
    x: u32,
    y: u32,
}

impl Cell {
    pub fn get_pos(&self) -> (u32, u32) {
        (self.x, self.y)
    }

    pub fn dist(&self, other: CellRef) -> f64 {
        let x_len = (self.x as i32 - other.x as i32).abs();
        let y_len = (self.y as i32 - other.y as i32).abs();
        ((x_len.pow(2) + y_len.pow(2)) as f64).sqrt()
    }
}

impl fmt::Display for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Customize so only `x` and `y` are denoted.
        write!(f, "({}, {})", self.x, self.y)
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn same_slope() {
        let canvas = Canvas::new();

        let true_cell = Rc::new(Cell { x: 2, y: 2 });
        let true_compare_cells = vec![Rc::new(Cell { x: 0, y: 0 }), Rc::new(Cell { x: 4, y: 4 })];
        // println!("{}", canvas.filter_same_slope(true_cell, true_compare_cells).is_some());
        // assert!(canvas.filter_same_slope(true_cell, true_compare_cells).is_some());

        let false_cell = Rc::new(Cell { x: 0, y: 0 });
        let false_compare_cells = vec![Rc::new(Cell { x: 4, y: 0 }), Rc::new(Cell { x: 0, y: 4 }), Rc::new(Cell { x: 4, y: 4 })];
        // assert!(canvas.filter_same_slope(false_cell, false_compare_cells).is_none());
    }
}
// #[derive(Eq, Clone)]
// pub struct Dist {
//     start: CellRef,
//     end: CellRef,
// }

// impl Dist {
//     pub fn get_legs(&self) -> (u32, u32) {
//         let x_leg = (self.start.x as i32 - self.end.x as i32).abs();
//         let y_leg = (self.start.y as i32 - self.end.y as i32).abs();
//         (x_leg as u32, y_leg as u32)
//     }

//     pub fn get_dists(&self) -> (u32, u32) {
//         let x_dist = self.start.x - self.end.x;
//         let y_dist = self.start.y - self.end.y;
//         (x_dist, y_dist)
//     }

//     pub fn get_mag(&self) -> f64 {
//         let (x_leg, y_leg) = self.get_legs();
//         ((x_leg.pow(2) + y_leg.pow(2)) as f64).sqrt()
//     }

//     pub fn cmp_slope(&self, other: &Self) -> bool {
//         let (x_leg_a, y_leg_a) = self.get_legs();

//         for x_offset in -1..2 {
//             for y_offset in -1..2 {
//                 let (x_leg_b, y_leg_b) = (other.get_legs().0 as i32 + x_offset, other.get_legs().1 as i32 + y_offset);


//                 if x_leg_a == 0 && x_leg_b == 0 {
//                     return true;
//                 } else if x_leg_a == 0 || x_leg_b == 0 {
//                     continue;
//                 }

//                 // 1st grade math :p
//                 let slope_a = y_leg_a as f64 / x_leg_a as f64;
//                 let slope_b = y_leg_b as f64 / x_leg_b as f64;

//                 if slope_a == slope_b {
//                     return true;
//                 }
//             }
//         }

//         false
//     }
// }

// impl Ord for Dist {
//     fn cmp(&self, other: &Self) -> Ordering {
//         // self.get_mag().cmp(&other.get_mag())
//         let (a, b) = (self.get_mag(), other.get_mag());
//         if a > b { Ordering::Greater }
//         else if a < b { Ordering::Less }
//         else { Ordering::Equal }
//     }
// }

// impl PartialOrd for Dist {
//     fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
//         Some(self.cmp(other))
//     }
// }

// impl PartialEq for Dist {
//     fn eq(&self, other: &Self) -> bool {
//         self.get_mag() == other.get_mag()
//     }
// }




