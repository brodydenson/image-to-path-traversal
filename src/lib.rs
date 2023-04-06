
// Storing Characters
// X: Timeline
// Y: X coordinates of character
// Z: Y coordinates of character
//
// ......
// ......
// ......
// ......
// ......
//
// Input:
// ..##@.
// .#....
// ..##..
// ....#.
// .###..
//
// ...1x0
// ...000
// ......
// ......
// ......
//
// ...x1.
// ......
// ......
// ......
// ......
//

use std::fmt;
use std::cmp::Ordering;
use std::rc::Rc;
use image::{ImageBuffer, Luma};
use std::path::Path;
type CellRef = Rc<Cell>;

pub struct Canvas {
    grid: Grid,
    outline: Vec<CellRef>,
    skeleton: Vec<CellRef>,
}


impl Canvas {
    pub fn new() -> Canvas {
        Canvas { grid: Grid::new(), outline: Vec::new(), skeleton: Vec::new() }
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

        println!("Building outline...");
        self.build_outline();
        println!("Building skeleton... (may take a while)");
        self.build_skeleton();
    }

    pub fn save(&self) {
        let original_path = Path::new("original.png");
        let outline_path = Path::new("outline.png");
        let skeleton_path = Path::new("skeleton.png");

        let (width, height) = self.grid.dim;
        // Saving original image
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            *pixel = match self.grid.at(x, y) {
                Some(_) => Luma([0]),
                _ => Luma([255]),
            }
        }

        new_img
            .save(original_path)
            .expect("Error saving original image");

        // Saving outline image
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            let in_outline = self.outline
                .iter()
                .find(|cell| cell.x == x as i32 && cell.y == y as i32);
            *pixel = match in_outline {
                Some(_) => Luma([0]),
                _ => Luma([255]),
            }
        }

        new_img
            .save(outline_path)
            .expect("Error saving outline image");

        // Saving skeleton image
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            let in_skeleton = self.skeleton
                .iter()
                .position(|cell| cell.x == x as i32 && cell.y == y as i32);
            let in_outline = self.outline
                .iter()
                .find(|cell| cell.x == x as i32 && cell.y == y as i32);
            *pixel = match in_skeleton {
                Some(i) => Luma([(194f64 * (i as f64 / self.skeleton.len() as f64)) as u8]),
                _ => match in_outline {
                    Some(_) => Luma([128]),
                    _ => Luma([225]),
                },
            }
        }

        new_img
            .save(skeleton_path)
            .expect("Error saving skeleton image");
    }

    fn build_outline(&mut self) {
        self.outline = self.grid.select(|cell| self.grid.get_adjacent(&cell).len() < 4);
    }

    fn build_skeleton(&mut self) {
        let mut skeleton: Vec<CellRef> = Vec::new();

        // Gotta silence the compiler somehow
        let outline = self.outline.clone();
        let interior = self.grid.select(|cell| !outline.contains(&Rc::clone(cell)));
        // let interior: Vec<CellRef> = taken_cells
        //     .into_iter()
        //     .filter(|cell| !outline.contains(&Rc::clone(cell)))
        //     .collect();

        // Creating interior cell dists
        let mut in_cell_dists: Vec<(CellRef, Dist)> = Vec::new();

        // Iterating over all interiors to compare with the outline
        for in_cell in interior.iter() {
            // Get all dists from in_cell with every outline
            let dists: Vec<Dist> = outline
                .iter()
                .map(|out_cell| Dist { start: Rc::clone(&in_cell), end: Rc::clone(out_cell) })
                .collect();

            let min_dist = dists.iter().min().unwrap();

            // Sorry rustaceans, but I have to use clone
            in_cell_dists.push((Rc::clone(&in_cell), min_dist.clone()));
        }

        in_cell_dists.sort_by(|(_, dist_a), (_, dist_b)| dist_b.cmp(dist_a));

        let (mut count, limit) = (0usize, 600000usize);
        while in_cell_dists.len() > 0 && count < limit {
            // if count == 590 {
            //     skeleton = vec![];
            // }
            let in_cell = &Rc::clone(&in_cell_dists[0].0);
            // println!("Interior Length: {}", in_cell_dists.len());

            // Get all dists from in_cell with every outline (again)
            // Optimization Idea: Don't do this...
            let mut dists: Vec<(CellRef, Dist)> = outline
                .iter()
                .map(|out_cell| (Rc::clone(out_cell), Dist{ start: Rc::clone(&in_cell), end: Rc::clone(out_cell) }))
                .collect();

            dists.sort_by(|(_, dist_a), (_, dist_b)| dist_a.cmp(dist_b));
            
            let find_sorrounding = |out_cell_a: &CellRef, out_dist_a: &Dist, targets: Vec<(i32, i32)>| {
                for (target_x, target_y) in targets {
                    let result = dists.iter().find(|(out_cell, _)| out_cell.get_pos() == (target_x, target_y));

                    if let Some((out_cell_b, out_dist_b)) = result {
                        let hypotenuse = Dist { start: Rc::clone(out_cell_a), end: Rc::clone(out_cell_b) };

                        // Making sure that the two outline points aren't right next to each other
                        if hypotenuse.get_mag() >= out_dist_a.get_mag() && hypotenuse.get_mag() >= out_dist_b.get_mag() {
                            return result;
                        }
                    }
                }
                None
            };

            let (min_out_cell_a, min_out_dist_a) = &dists[0];
            let (x_dist, y_dist) = min_out_dist_a.get_dists();
            'outer: for x_offset in -1..2 {
                for y_offset in -1..2 {
                    let targets = vec![
                        (in_cell.get_pos().0 + (x_dist + x_offset), in_cell.get_pos().1 + (y_dist + y_offset)),
                        // (in_cell.get_pos().0 + (-x_dist + x_offset), in_cell.get_pos().1 + (y_dist + y_offset)),
                        // (in_cell.get_pos().0 + (x_dist + x_offset), in_cell.get_pos().1 + (-y_dist + y_offset)),
                    ];

                    let result = find_sorrounding(min_out_cell_a, min_out_dist_a, targets);

                    if let Some((min_out_cell_b, min_out_dist_b)) = result {
                        println!("Cell A: {}, Cell B (target): {}", min_out_cell_a, min_out_cell_b);


                        // println!("Before length: {}", in_cell_dists.len());
                        self.erase_in_cells(&mut in_cell_dists, Rc::clone(in_cell), Rc::clone(&min_out_cell_a));
                        self.erase_in_cells(&mut in_cell_dists, Rc::clone(in_cell), Rc::clone(&min_out_cell_b));
                        // println!("After length: {}", in_cell_dists.len());

                        skeleton.push(Rc::clone(&in_cell));
                        skeleton.push(Rc::clone(&min_out_cell_a));
                        skeleton.push(Rc::clone(&min_out_cell_b));
                        break 'outer;
                    }
                }
            }
            

            // let min_leg_sum = &dists[0].1.get_legs().0 + &dists[0].1.get_legs().1;
            //
            //
            // println!("\nNew sum: {}", min_leg_sum);
            // while dists.len() >= 2 {
            //     let pair = (&dists[0], &dists[1]);
            //     let (legs_a, legs_b) = (pair.0.1.get_legs(), pair.1.1.get_legs());
            //     let (out_cell_a, out_cell_b) = (&pair.0.0, &pair.1.0);
            //_b
            //     println!("Positions: ");
            //     println!("{}, {}", out_cell_a.get_pos().0, out_cell_a.get_pos().1);
            //     println!("{}, {}", out_cell_b.get_pos().0, out_cell_b.get_pos().1);
            //
            //     // if legs_a.0 as i32 - 1 > legs_a.0 as i32 && legs_a.1 as i32 - 1 > legs_b.1 as i32 {
            //     //     break;
            //     // }
            //
            //     // Letting room for a 1px margin of error
            //     if (legs_a.0 + legs_a.1) as i32 - 2 > min_leg_sum as i32 || (legs_b.0 + legs_b.1) as i32 - 2 > min_leg_sum as i32 {
            //         // println!("{}: ({}, {}), ({}, {})", min_leg_sum, legs_a.0, legs_a.1, legs_b.0, legs_b.1);
            //         println!("Too far");
            //         dists.remove(1);
            //         continue;
            //     }
            //
            //     if !pair.0.1.cmp_slope(&pair.1.1) {
            //         println!("Different slope");
            //         // dists.remove(0);
            //         dists.remove(1);
            //         continue;
            //     }
            //
            //     let hypotenuse = Dist { start: Rc::clone(out_cell_a), end: Rc::clone(out_cell_b) };
            //
            //     // Making sure that the two outline points aren't right next to each other
            //     if hypotenuse.get_mag() < pair.0.1.get_mag() || hypotenuse.get_mag() < pair.1.1.get_mag() {
            //         println!("Too close");
            //         dists.remove(1);
            //         // dists.remove(0);
            //         continue;
            //     }
            //
            //     println!("Passed");
            //
            //
            //     // println!("Before length: {}", in_cell_dists.len());
            //     self.erase_in_cells(&mut in_cell_dists, Rc::clone(in_cell), Rc::clone(&out_cell_a));
            //     self.erase_in_cells(&mut in_cell_dists, Rc::clone(in_cell), Rc::clone(&out_cell_b));
            //     // println!("After length: {}", in_cell_dists.len());
            //
            //     skeleton.push(Rc::clone(&in_cell));
            //     skeleton.push(Rc::clone(&out_cell_a));
            //     skeleton.push(Rc::clone(&out_cell_b));
            //     break;
            // }


            // Find and remove element to not be reused
            if let Some(index) = in_cell_dists.iter().position(|(cell, _)| *cell == Rc::clone(in_cell)) {
                // Remove interior and in_cell_dists because in_ceLl_dists is being accessed
                // inside this loop
                // interior.remove(index);
                in_cell_dists.remove(index);
            }

            // in_cell_dists.remove(0);

            // println!("Iteration Done.");
            count += 1;
        }

        // for i in &skeleton {
        //     println!("{}, {}", i.get_pos().0, i.get_pos().1);
        // }

        self.skeleton = skeleton;
    }

    fn connect_corners(&mut self, skeleton: &mut Vec<CellRef>) {

    }


    // Brensenham's Line algorithm
    fn erase_in_cells(&mut self, in_cell_dists: &mut Vec<(CellRef, Dist)>, c1: CellRef, c2: CellRef) {
        let mut remove_in_cell = |x, y| {
            // Biggest hault
            let half_thickness = (1.0 / 2.0) as isize;
            for dx in -half_thickness..=half_thickness {
                for dy in -half_thickness..=half_thickness {
                    if dx * dx + dy * dy <= half_thickness * half_thickness {
                        if let Some(index) = in_cell_dists.iter().position(|(cell, _)| cell.get_pos() == (x + dx as i32, y + dy as i32)) {
                            // self.skeleton.push(Rc::clone(&in_cell_dists[index].0));
                            in_cell_dists.remove(index);
                        }
                    }
                }
            }
        };

        let (x1, y1) = c1.get_pos();
        let (x2, y2) = c2.get_pos();

        let dx = (x2 - x1).abs();
        let dy = (y2 - y1).abs();
        let (mut x, mut y) = (x1, y1);
        let sx = if x1 < x2 { 1 } else { -1 };
        let sy = if y1 < y2 { 1 } else { -1 };

        if dx > dy {
            let mut err = dx / 2;
            while x != x2 {
                remove_in_cell(x, y);
                err -= dy;
                if err < 0 {
                    y += sy;
                    err += dx;
                }
                x += sx;
            }
        } else {
            let mut err = dy / 2;
            while y != y2 {
                remove_in_cell(x, y);
                // points.push((x, y));
                err -= dx;
                if err < 0 {
                    x += sx;
                    err += dy;
                }
                y += sy;
            }
        }

        remove_in_cell(x, y);
    }

    pub fn print(&self) {
        self.grid.print();
    }

    pub fn print_outline(&self) {
        self.grid.print_vertices(&self.outline);
    }

    pub fn print_skeleton(&self) {
        self.grid.print_vertices(&self.skeleton);
    }
}

#[derive(Eq, Clone)]
pub struct Dist {
    start: CellRef,
    end: CellRef,
}

impl Dist {
    pub fn get_legs(&self) -> (u32, u32) {
        let x_leg = (self.start.x - self.end.x).abs();
        let y_leg = (self.start.y - self.end.y).abs();
        (x_leg as u32, y_leg as u32)
    }

    pub fn get_dists(&self) -> (i32, i32) {
        let x_dist = self.start.x - self.end.x;
        let y_dist = self.start.y - self.end.y;
        (x_dist, y_dist)
    }

    pub fn get_mag(&self) -> f64 {
        let (x_leg, y_leg) = self.get_legs();
        ((x_leg.pow(2) + y_leg.pow(2)) as f64).sqrt()
    }

    pub fn cmp_slope(&self, other: &Self) -> bool {
        // Making sure that the points aren't right next to each other
        // let end_dist = Dist { start: Rc::clone(&self.end), end: Rc::clone(&other.end) };
        // if end_dist.get_mag() <= 2f64.sqrt() {
        //     println!("Too close: (Start: ({}, {}), End: ({}, {})) (Start: ({}, {}), End: ({}, {})",
        //         self.start.get_pos().0, self.start.get_pos().1, self.end.get_pos().0, self.end.get_pos().1,
        //         other.start.get_pos().0, other.start.get_pos().1, other.end.get_pos().0, other.end.get_pos().1);
        //     return false;
        // }

        let (x_leg_a, y_leg_a) = self.get_legs();

        for x_offset in -1..2 {
            for y_offset in -1..2 {
                let (x_leg_b, y_leg_b) = (other.get_legs().0 as i32 + x_offset, other.get_legs().1 as i32 + y_offset);


                if x_leg_a == 0 && x_leg_b == 0 {
                    return true;
                } else if x_leg_a == 0 || x_leg_b == 0 {
                    continue;
                }

                // 1st grade math :p
                let slope_a = y_leg_a as f64 / x_leg_a as f64;
                let slope_b = y_leg_b as f64 / x_leg_b as f64;

                if slope_a == slope_b {
                    return true;
                }
            }
        }

        false
    }
}

impl Ord for Dist {
    fn cmp(&self, other: &Self) -> Ordering {
        // self.get_mag().cmp(&other.get_mag())
        let (a, b) = (self.get_mag(), other.get_mag());
        if a > b { Ordering::Greater }
        else if a < b { Ordering::Less }
        else { Ordering::Equal }
    }
}

impl PartialOrd for Dist {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Dist {
    fn eq(&self, other: &Self) -> bool {
        self.get_mag() == other.get_mag()
    }
}


pub struct Grid {
    dim: (u32, u32),
    cells: Vec<Vec<Option<CellRef>>>,
}
 
impl Grid {
    pub fn new() -> Grid {
        Grid { dim: (0, 0), cells: Vec::new() }
    }

    pub fn set_dim(&mut self, dim: (u32, u32)) {
        self.dim = dim;
        self.cells = vec![vec![None; dim.1 as usize]; dim.0 as usize];
    }

    pub fn set_grid(&mut self, str_grid: &str) {
        let str_grid_fmt: String = str_grid.chars().filter(|c| !c.is_whitespace()).collect();
        let (height, width) = self.dim;

        for row in 0..height {
            let mut row_vec: Vec<Option<CellRef>> = Vec::new();
            for col in 0..width {
                let cell = match str_grid_fmt.chars().nth((row * width + col) as usize).unwrap() {
                    // '#' => Some(Rc::new(RefCell::new(Cell { x: row as i32, y: col as i32 }))),
                    '#' => Some(Rc::new(Cell { x: row as i32, y: col as i32 })),
                    '.' => None,
                    c => panic!("Invalid symbol {}", c),
                };
                row_vec.push(cell);
            }
            self.cells.push(row_vec);
        }
    }

    pub fn set_cell(&mut self, row: u32, col: u32, some: bool) {
        let cell_opt = if some {
            Some(Rc::new(Cell { x: row as i32, y: col as i32 }))
        } else {
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

        let find = |x: i32, y: i32| {
            for v in vertices {
                if v.get_pos() == (x, y) {
                    return true;
                }
            }
            false
        };

        for row in 0..height {
            for col in 0..width {
                if find(row as i32, col as i32) {
                    print!("#");
                } else {
                    print!(".");
                }
            }
            print!("\n");
        }
    }

    pub fn at(&self, row: u32, col: u32) -> Option<CellRef> {
        // if let Some(row_vec) = self.cells.get(row) {
        //     if let Some(cell_opt) = row_vec.get(col) {
        //         // return Some(cell_opt.as_ref().map(Rc::clone));
        //         return Rc::new(cell_opt);
        //     }
        // }
        // None
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
                // Adjecent not including the element itself, and diaginal
                if (row_offset != 0 && col_offset != 0)
                    || (row_offset == 0 && col_offset == 0) {
                    continue;
                }

                let row = x + row_offset;
                let col = y + col_offset;

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
    x: i32,
    y: i32,
}

impl Cell {
    pub fn get_pos(&self) -> (i32, i32) {
        (self.x, self.y)
    }

    pub fn dist(&self, other: CellRef) -> f64 {
        let x_len = (self.x - other.x).abs();
        let y_len = (self.y - other.y).abs();
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
