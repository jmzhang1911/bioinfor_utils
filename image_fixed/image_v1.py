import sys

sys.path.append(r'D:\workspace\MyUtils')
from myrunner import MyPath
from pathlib import Path
from tqdm import tqdm
import numpy as np
import cv2 as cv
import argparse
import logging
import math


class MyImageTools:
    @staticmethod
    def cat_pic(image_obj):
        dim = image_obj.size
        v_concat_list = []
        for line in range(dim[0]):
            h_concat_list = []
            for pos, piece in image_obj.container.items():
                if pos[1] == line:
                    h_concat_list.append(piece)
                h_concat_list = sorted(h_concat_list, key=lambda x: x.relative_position[0])
                v_concat_list.append(cv.hconcat(h_concat_list))
        real_plot = cv.vconcat([v_concat_list])
        cv.imwrite('real_plot.tiff', real_plot)

    @staticmethod
    def order_dict(raw_dict: dict, dim: tuple):
        start = 0
        x_pos = 0
        all_position = [i for i in raw_dict.keys()]
        ordered_dict = {}

        for x_lim in range(dim[0]):
            b = sorted(all_position, key=lambda x: x[0])[start:start + dim[1]]
            c = sorted(b, key=lambda x: x[1])
            y_pos = 0
            for i in c:
                ordered_dict[(x_pos, y_pos)] = raw_dict[i]
                y_pos += 1

            start += dim[1]
            x_pos += 1

        return ordered_dict

    @staticmethod
    def adjust_piece(*piece_path, dim=(46, 46)):
        # for piece in piece_path:
        # MyImage中同一个位置的MyPiece，取这些位置中score为2的，若没有score为2的取，点数最接近1050的
        MyPath.mkdir('right_pieces')
        right_image = MyImage()
        pos = [(x, y) for x in range(dim[0]) for y in range(dim[1])]
        piece_in_all = {}
        for x, y in tqdm(pos, desc='picking the right pieces'):
            for piece in piece_path:
                piece_plot = piece + '/piece_{}_{}.tiff'.format(str(x), str(y))
                if list(Path(piece_plot).stat())[6] < 400000:
                    continue
                plot = cv.imread(str(piece_plot))
                piece_obj = MyPiece(raw_plot=plot, relative_position=(x, y))
                piece_in_all[Path(piece).stem] = piece_obj

            score_list = [piece_obj for file_name, piece_obj in piece_in_all.items() if piece_obj.score == 2]
            if score_list:
                right_image.container[(x, y)] = score_list[0]
                # cv.imwrite('right_pieces/piece_{}_{}.tiff'.format(str(x), str(y)), score_list[0].fixed_plot)
            else:
                spot_num = [(piece_obj, math.fabs(piece_obj.con_number - 1050)) for file_name, piece_obj in
                            piece_in_all.items()]
                closest = sorted(spot_num, key=lambda x: x[1])[0][0]
                right_image.container[(x, y)] = closest
                # cv.imwrite('right_pieces/piece_{}_{}.tiff'.format(str(x), str(y)), closest.fixed_plot)

        MyImageTools.cat_pic(right_image)


class MySpot:
    def __init__(self, center_position: tuple, relative_position=None, msg=None, con=None):
        self.center_position = center_position  # 中心点的绝对坐标
        self.relative_position = relative_position  # spot在spices中的相对位置，即排序之后(x, y)
        self.containers = [msg]  # 点的信息，即荧光信号
        self.con = con  # 点的边框（相较于MyPiece对象而言），array


class MyPiece:
    def __init__(self, raw_plot, center_position=None, relative_position=None, size=(30, 35)):
        self.center_position = center_position  # 中心点在完整图片中的绝对坐标
        self.relative_position = relative_position  # spot在spices中的相对位置，即排序之后(x, y)
        self.raw_plot = raw_plot  # piece的原图片
        self.shape = raw_plot.shape  # piece原图片的尺寸466 * 467 (行，列)
        self.size = size  # piece行列点数

        self.fixed_plot = None  # 校正后的图
        self.containers = {}  # piece容器，存放点的信息，即{(x, y): MySpotObject ...}
        self.score = None  # 此piece的完整得分，0:不完整，1:基本完整，2:非常完整
        self.score_precise = None  # 精确打分
        self.con = None  # 存储所有spot边框信息
        self.con_number = None

        # 添加属性：解析piece中的点，调整raw_plot的尺寸得到fixed_plot，根据点的信息对piece进行打分
        self.resolve_to_spot()
        # self.score_it_precise()

    @staticmethod
    def read_piece(piece, center_position=None, relative_position=None, size=(30, 35)):
        p = cv.imread(str(piece))
        piece_obj = MyPiece(raw_plot=p, center_position=center_position,
                            relative_position=relative_position, size=size)

        return piece_obj

    def __add__(self, other):
        """两个多个pieces相加，将pieces中相同位置的点的信息相加，返回一个新的Piece"""
        return self

    @staticmethod
    def con_fixed(raw_con):
        """修复两两相连点，将其拆分为两个"""
        new_con = []
        for con in raw_con:
            x_list = [x for [[x, _]] in con]
            y_list = [y for [[_, y]] in con]
            x_res = [max(x_list), min(x_list)]
            y_res = [max(y_list), min(y_list)]

            x_len = x_res[0] - x_res[1]
            y_len = y_res[0] - y_res[1]

            # 将连在一起的点分开
            if 80 <= cv.contourArea(con) <= 200:
                # 横向：
                if x_len > y_len:
                    con_1 = np.array([x >= (x_res[1] + 4) and [[x_res[1] + 8, y]] or [[x, y]] for [[x, y]] in con])
                    con_2 = np.array([x <= (x_res[0] - 4) and [[x_res[0] - 8, y]] or [[x, y]] for [[x, y]] in con])

                else:
                    con_1 = np.array(
                        [y >= (y_res[1] + 4) and [[x_res[1] + 8, y_res[1] + 4]] or [[x, y]] for [[x, y]] in con])
                    con_2 = np.array(
                        [y <= (y_res[0] - 4) and [[x_res[0] - 8, y_res[0] - 4]] or [[x, y]] for [[x, y]] in con])

                new_con.append(con_1)
                new_con.append(con_2)
            else:
                new_con.append(con)

        return new_con

    def resolve_to_spot(self):
        """解析piece中的spot，返回一个字典存储{(x, y): MySpotObject ...}暂时存放在list中"""
        p = self.raw_plot.copy()
        gray = cv.cvtColor(p, cv.COLOR_BGR2GRAY)  # 转换为灰度
        r, p3 = cv.threshold(gray, 200, 250, cv.THRESH_BINARY_INV)  # 颜色二值化
        p3 = cv.medianBlur(p3, 7)  # 模糊
        binary, contours, hierarchy = cv.findContours(p3, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_NONE)
        contours = [con for con in contours if cv.contourArea(con) >= 5]

        new_con = self.con_fixed(contours)
        new_con = [con for con in new_con if cv.contourArea(con) <= 120]

        self.con = new_con  # MyPiece的con
        self.adjust_it()  # 基于调整后的con进行拆分
        self.con_number = len(self.con)  # 添加此piece的个数，

        piece_dict = {}
        for con in new_con:
            x_list = [x for [[x, _]] in con]
            y_list = [y for [[_, y]] in con]
            x_res = [max(x_list), min(x_list)]
            y_res = [max(y_list), min(y_list)]
            center_x = x_res[1] + round((x_res[0] - x_res[1]) * 0.5, 3)
            center_y = y_res[1] + round((y_res[0] - y_res[1]) * 0.5, 3)

            # 点的信息
            spot_msg = ''
            self.score_it()

            if self.score == 2:
                piece_dict[(center_x, center_y)] = ((center_x, center_y), spot_msg, con)

        # 如果点的个数为1050个则构建MySpot对象
        piece_dict = MyImageTools.order_dict(piece_dict, dim=(30, 35))
        for pos, value in piece_dict.items():
            spot_obj = MySpot(center_position=value[0], msg=value[1], relative_position=pos, con=value[2])
            self.containers[pos] = spot_obj

    def adjust_it(self, gap=470):
        """
        裁剪，统一所有的大小，方便piece之间的拼接：思路确定所有con中四个方向的最大值，结果保存于fixed_plot中
        填充比例为(480,480)
        """
        h, w, _ = self.raw_plot.shape
        new_con = [con for con in self.con if cv.contourArea(con) >= 10]  # 取面积大于10的con来确定piece周边
        x_res = []
        y_res = []
        for con in new_con:
            x_list = [x for [[x, _]] in con]
            y_list = [y for [[_, y]] in con]
            x_res += [max(x_list), min(x_list)]
            y_res += [max(y_list), min(y_list)]

        # logging.info('原始长宽{},{}'.format(h, w))
        # logging.info('点的下上左右{},{},{},{}'.format(min(y_res), max(y_res), min(x_res), max(x_res)))

        # 按照上下左右的最大值分别延申4个像素进行裁剪
        x_min = min(x_res) - 4 if min(x_res) - 4 >= 0 else 0  # 确定左边裁剪位置
        x_max = max(x_res) + 4 if max(x_res) + 4 <= w else w  # 确定右边裁剪位置
        y_min = min(y_res) - 4 if min(y_res) - 4 >= 0 else 0  # 确定下边裁剪位置
        y_max = max(y_res) + 4 if max(y_res) + 4 <= h else h  # 确定上边裁剪位置
        x_len, y_len = x_max - x_min, y_max - y_min  # 裁切矩形的长于狂，(width, height)

        # logging.info('移动后的下上左右{},{},{},{}'.format(y_min, y_max, x_min, x_max))

        # 与 470 * 470标准矩形缝隙，在不存在交集的情况下可以直接加减(可取负数)
        up, down = math.ceil((gap - y_len) * 0.5), math.floor((gap - y_len) * 0.5)  # 上下
        left, right = math.ceil((gap - x_len) * 0.5), math.floor((gap - x_len) * 0.5)  # 左右

        # logging.info('上下差：{}{}'.format(up, down))
        # logging.info('点的长宽{}:{}'.format(x_len, y_len))
        # logging.info('白边上下左右{},{},{},{}'.format(up, down, left, right))

        x_min_scale, x_max_scale = x_min - left, x_max + right  # x左右移动
        y_min_scale, y_max_scale = y_min - down, y_max + up  # y上下移动

        if (x_len <= gap and y_len <= gap) or (x_len > gap and y_len > gap):
            """正常缩放，完全包括或者被包括于480*480"""
            # logging.info('1')
            if x_len <= gap:
                # logging.info('a')
                logging.info('{},{},{},{}'.format(y_min, y_max, x_min, x_max))
                scale_plot = cv.copyMakeBorder(self.raw_plot[y_min:y_max, x_min:x_max],
                                               up, down, left, right, cv.BORDER_CONSTANT,
                                               value=(255, 255, 255))
            else:
                # logging.info('b')
                scale_plot = self.raw_plot[y_min_scale:y_max_scale, x_min_scale:x_max_scale]

        elif x_len > gap:
            # logging.info('2')
            """首先上下填充为完全包括480*480，左右裁剪成480"""
            up_spr, down_spr = math.ceil((gap - y_len) * 0.5), math.floor((gap - y_len) * 0.5)
            scale_plot = cv.copyMakeBorder(self.raw_plot[y_min:y_max, x_min_scale:x_max_scale],
                                           up_spr, down_spr, 0, 0, cv.BORDER_CONSTANT,
                                           value=(255, 255, 255))
            y_min_scale, y_max_scale = 0, gap
            # scale_plot = scale_plot[y_min_scale:y_max_scale, x_min_scale:x_max_scale]

        else:
            # logging.info('3')
            """首先左右填充为完全包括480*480，左右裁剪成480"""
            left_spr, right_spr = math.ceil((gap - x_len) * 0.5), math.floor((gap - x_len) * 0.5)

            # logging.info('左右阔{},{}'.format(left_spr, right_spr))

            scale_plot = cv.copyMakeBorder(self.raw_plot[y_min_scale:y_max_scale, x_min:x_max],
                                           0, 0, left_spr, right_spr, cv.BORDER_CONSTANT,
                                           value=(255, 255, 255))
            x_min_scale, x_max_scale = 0, gap

        # logging.info('上下裁{}, {}'.format(y_min_scale, y_max_scale))
        # scale_plot = scale_plot[y_min_scale:y_max_scale, x_min_scale:x_max_scale]

        self.fixed_plot = scale_plot

    def fix_it(self):
        """修复，待实现"""
        return self

    def score_it(self):
        """仅基于piece中点的个数对piece的完整度进行打分"""
        if self.con_number == 1050:
            self.score = 2
        elif math.fabs(self.con_number - 1050) <= 10:
            self.score = 1
        else:
            self.score = 0

    def score_it_precise(self):
        """
        - 解决点个数满足1050个但仍然是缺损的情况，暂时基于两个特征来判断
            - 基于横纵方向点的平均面积判断：正常情况面积的均值为估摸20~40
            - 基于横纵方向相邻点的距离判断
        """
        if self.score == 2:
            # 横向取所有点对象，(x, y)，一排30点
            for y in range(35):  # 可以取出35排
                # self.containers: {(x, y):MySpotObj,...}
                spot_list = [spot_obj for _, spot_obj in self.containers.items() if _[1] == y]
                # 对此横向的所有点排序，按照相对路径从小到大，用于计算点之间距离
                spot_list = sorted(spot_list, key=lambda i: i.relative_position[0])

                if len(spot_list) != 30:
                    logging.info('{}-{} wrong length y_len={} -> {}'.
                                 format(*self.relative_position, y, len(spot_list)))
                    self.score_precise = None
                    return

                # 计算此横横向所有点的平均面积
                spot_mean_area = np.mean([cv.contourArea(spot.con) for spot in spot_list])

                if spot_mean_area < 20 or spot_mean_area > 40:
                    logging.info('{}-{} wrong area in y={} -> {}'.
                                 format(*self.relative_position, y, spot_mean_area))
                    self.score_precise = 'bad'
                    return

                start, end = 0, 2
                for time in range(30 - 1):
                    spot_start, spot_end = spot_list[start:end]
                    (start_x, start_y) = spot_start.center_position
                    (end_x, end_y) = spot_end.center_position
                    distance = round(math.sqrt(math.pow((end_x - start_x), 2) + math.pow((end_y - start_y), 2)), 2)

                    # 经验值，后续需要调整
                    if distance < 10 or distance > 25:
                        logging.info('{}-{} wrong distance y:{}={}-{} -> {}'.
                                     format(*self.relative_position, y, start, end, distance))
                        self.score_precise = 'bad'
                        return

                    start += 1
                    end += 1

            # 纵向取所有点对象，(x, y)，一排35点，纵向对象存在交叉，及计数列和偶数列
            for x in range(30):  # 可以取出30排
                # self.containers: {(x, y):MySpotObj,...}
                spot_list = [spot_obj for _, spot_obj in self.containers.items() if _[0] == x]
                # 对此横向的所有点排序，按照相对路径从小到大，用于计算点之间距离
                spot_list = sorted(spot_list, key=lambda i: i.relative_position[1])
                spot_list_even = spot_list[1::2]
                spot_list_odd = spot_list[::2]

                if len(spot_list) != 35:
                    logging.info('{}-{} wrong length x={} -> {}'.
                                 format(*self.relative_position, x, len(spot_list)))
                    self.score_precise = None
                    return

                # 对此
                spot_mean_area_even = np.mean([cv.contourArea(spot.con) for spot in spot_list_even])
                spot_mean_area_odd = np.mean([cv.contourArea(spot.con) for spot in spot_list_odd])

                # logging.info('x={} odd={},even={}'.format(x, spot_mean_area_odd, spot_mean_area_even))

                if spot_mean_area_even < 20 or spot_mean_area_even > 40 or \
                        spot_mean_area_odd < 20 or spot_mean_area_odd > 40:
                    wrong_type = 'mean' if spot_mean_area_even < 20 else 'odd'
                    value = spot_mean_area_even if spot_mean_area_even < 20 else spot_mean_area_odd
                    logging.info('{}-{} wrong area type={} x={} -> {}'.
                                 format(*self.relative_position, wrong_type, x, value))
                    self.score_precise = 'bad'
                    return

                start, end = 0, 2
                for time in range(35 - 1):
                    spot_start, spot_end = spot_list[start:end]
                    (start_x, start_y) = spot_start.center_position
                    (end_x, end_y) = spot_end.center_position

                    distance = math.sqrt(math.pow((end_x - start_x), 2) + math.pow((end_y - start_y), 2))
                    distance = round(distance, 2)

                    # 经验值，后续需要调整
                    if distance < 13 or distance > 18:
                        logging.info('{}-{} wrong distance x:{}={}-{} -> {}'.
                                     format(*self.relative_position, x, start, end, distance))
                        self.score_precise = 'bad'
                        return

                start += 1
                end += 1

            self.score_precise = 'good'

    def draw_con(self, pos=None, plot=None):
        """
        因为是在raw_plot的基础上进行的con识别，因此将con画在raw_plot上，pos为一个tuple(x, y)坐标，可叠加绘图传入tuple(plot,)
        (None,3):返回第三行
        (3,None):返回第三列
        (3,3):返回第三列第三个
        """
        p = self.raw_plot.copy() if not plot else plot[0].copy()
        if pos:
            x, y = pos
            if x and y:
                logging.info('1')
                con_list = [spot_obj.con for pos, spot_obj in self.containers.items() if pos[0] == x and pos[1] == y]
            elif x:
                con_list = [spot_obj.con for pos, spot_obj in self.containers.items() if pos[0] == x]
            elif y:
                con_list = [spot_obj.con for pos, spot_obj in self.containers.items() if pos[1] == y]
            else:
                con_list = self.con

            return cv.drawContours(p, con_list, -1, (0, 0, 255), 1)

        return cv.drawContours(p, self.con, -1, (0, 0, 255), 1)


class MyImage:
    kernel = np.ones((3, 3), np.uint8)  # 膨胀
    kernel2 = np.ones((10, 10), np.uint8)  # 腐蚀

    def __init__(self, plot_file=None, size=(46, 46)):
        if plot_file:
            self.plot_name = Path(plot_file).stem  # 图片名称，取stem
            self.plot = cv.imread(str(plot_file))  # 图片原图

            # 结果目录
            self.results_dir = Path(self.plot_name + '_image_results')  # 单张图片的结果目录
            MyPath.mkdir(str(self.results_dir))

        self.size = size  # 图片的行列piece数目(行 * 列)
        self.container = {}  # 存储piece {(x, y): MyPieceObject ...}

        # 边框矩形，拆分矩形
        self.piece_number = None  # 包含的piece个数
        self.con = []  # piece的边框
        self.rect = []  # 拆分用的矩形

    @staticmethod
    def _con_fixed(raw_con):
        """对矩识别到的矩形边框进行修正，此方法只可以修正横竖两两相连的情况"""
        rect_list = []
        for con in raw_con:
            x, y, w, h = cv.boundingRect(con)
            x, y, w, h = x, y, x + w, y + h
            x_len, y_len = w - x, h - y
            if x_len * y_len <= 60000:
                continue
            elif x_len > 2 * y_len:
                rect_list.append((w - y_len, y, w, h))
                rect_list.append((x, y, x + y_len, h))
            elif y_len > 2 * x_len:
                rect_list.append((x, y, w, y + x_len))
                rect_list.append((x, h - x_len, w, h))
            else:
                rect_list.append((x, y, w, h))

        return rect_list

    def get_con(self):
        p = self.plot.copy()
        gray = cv.cvtColor(p, cv.COLOR_BGR2GRAY)
        r, p3 = cv.threshold(gray, 250, 255, cv.THRESH_BINARY_INV)

        # 添加一次模糊origin中是没有的
        # p3 = cv.medianBlur(p3, 5)

        opening = cv.dilate(p3, self.kernel, iterations=2)  # origin = 2
        opening2 = cv.erode(opening, self.kernel2, iterations=10)  # origin = 10
        binary, contours, hierarchy = cv.findContours(opening2, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_NONE)

        # 获取腐蚀后边框的矩形边框
        new_cons = [cv.approxPolyDP(con, 0.1 * cv.arcLength(con, True), True) for con in contours]
        # 初步过滤一下面积较小的矩形边框
        new_cons = [con for con in new_cons if cv.contourArea(con) >= 10000]
        self.con = self._con_fixed(new_cons)
        # 统计得到的修正后的边框个数，目前个数为2116
        self.piece_number = str(len(self.con))
        logging.info('con number is ' + self.piece_number)

        return self.con

    def draw_con(self, save_plot_name='con_plot.tiff'):
        """在完整图形中绘制边框"""
        logging.info('drawing total con')
        p = self.plot.copy()
        for num, pos in enumerate(self.con):
            x, y, w, h = pos
            mean_x, mean_y = math.ceil(x + (w - x) / 2), math.ceil(y + (h - y) / 2),
            p = cv.rectangle(p, (x, y), (w, h), (0, 255, 0), 3)
            p = cv.putText(p, 'pos={}'.format(num),
                           (mean_x, mean_y), cv.FONT_HERSHEY_SIMPLEX, 1.2,
                           (255, 0, 0), 2)
        logging.info('saving con_plot.tiff ...')
        cv.imwrite(str(self.results_dir / save_plot_name), p)

    def draw_rect(self, save_plot_name='total_pieces.tiff'):
        """在完整大图中绘制拆分矩形"""
        logging.info('drawing total pieces...')
        p = self.plot.copy()
        for rect in self.rect:
            x, y, w, h = rect
            p = cv.rectangle(p, (x, y), (w, h), (0, 255, 0), 3)

        cv.imwrite(str(self.results_dir / save_plot_name), p)

    def split_by_con(self, pieces_dir='all_pieces', to_save=True):
        """
        实现拆分，构建拆分后的piece对象存入container中
        """
        p = self.plot.copy()
        plot_dict = {}
        for rect in self.con:
            x, y, w, h = rect
            x_left = math.ceil((w - x) * 0.17)  # 左点左移
            x_up = math.ceil((w - x) * 0.16)  # 左点上移
            y_right = math.ceil((h - y) * 0.13)  # 右点右移
            y_down = math.ceil((h - y) * 0.12)  # 右点下移

            x = x - x_left
            y = y - x_up
            w = w + y_right
            h = h + y_down

            # 截取保存
            center_x, center_y = math.ceil((w - x) * 0.5 + x), math.ceil((h - y) * 0.5 + h)
            plot_dict[(center_x, center_y)] = (p[y:h, x:w], (center_x, center_y))
            self.rect.append((x, y, w, h))

        # 截取，排序，保存
        res_path = str(self.results_dir / pieces_dir)
        plot_dict = MyImageTools.order_dict(plot_dict, dim=(46, 46))

        MyPath.mkdir(res_path)
        MyPath.mkdir(str(self.results_dir / 'good_pieces'))
        MyPath.mkdir(str(self.results_dir / 'ambiguous'))

        # 保存图片，将构建MySpot对象并存入MyImage的container中
        for k, plot in tqdm(plot_dict.items(), desc='saving and making pieces...'):
            p, (center_x, center_y) = plot
            piece_obj = MyPiece(center_position=(center_x, center_y), relative_position=k, raw_plot=p)
            self.container[k] = piece_obj
            if to_save:
                cv.imwrite('{}/piece_{}_{}.png'.format(res_path, *k), piece_obj.fixed_plot)

                # 2分的图片为good，保存2分的图片
                if piece_obj.score == 2:
                    cv.imwrite('{}/good_pieces/piece_{}_{}.png'.
                               format(self.results_dir, *k), piece_obj.fixed_plot)
                else:
                    cv.imwrite('{}/ambiguous/piece_{}_{}.png'.
                               format(self.results_dir, *k), piece_obj.fixed_plot)

    @staticmethod
    def read_pieces(pieces_dir):
        my_image = MyImage()
        pieces_file = [file for file in Path(pieces_dir).glob('*.png')]
        for piece in tqdm(pieces_file, desc='reading pieces'):
            p = cv.imread(str(piece))
            relative_x, relative_y = piece.stem.split('_')[1:3]
            relative_pos = (int(relative_x), int(relative_y))
            my_image.container[relative_pos] = MyPiece(raw_plot=p, relative_position=relative_pos)

        return my_image


if __name__ == '__main__':
    # 测试：
    pic_list = [r'D:\workspace\MyUtils\fun\v1\test\20211028-lcy-20X-1_Wholeslide_Default_Extended.tif',
                r'D:\workspace\MyUtils\fun\v1\test\20211028-lcy-20X-2_Wholeslide_Default_Extended.tif',
                r'D:\workspace\MyUtils\fun\v1\test\20211028-lcy-20X-3_Wholeslide_Default_Extended.tif',
                r'D:\workspace\MyUtils\fun\v1\test\20211028-lcy-20X-4_Wholeslide_Default_Extended.tif',
                r'D:\workspace\MyUtils\fun\v1\test\20211028-lcy-20X-5_Wholeslide_Default_Extended.tif']

    parser = argparse.ArgumentParser(description='image fixer', formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers()

    # 1)第一步：拆分图片得到piece存放于all_pieces，并将图片分为good_pieces和ambiguous_pieces
    image_splitter = subparsers.add_parser('step1_splitter', help='split image and divided into good/ambiguous')
    image_splitter.add_argument('-i', '--input', type=str, help='input images seperated by "," ps:image1,image2 ...)')
    image_splitter.set_defaults(action=('step1_splitter', lambda: None))

    # 2)第二步：读入所有图片的good_pieces将其拼装成完整图片
    image_concat = subparsers.add_parser('step2_concat', help='image concat')
    image_concat.add_argument('-g', '--good_pieces', type=str, help='all good pieces dir')
    image_concat.add_argument('-r', '--results', type=str, help='all good pieces dir', default='bingo')
    image_concat.set_defaults(action=('step2_concat', lambda: None))

    # 解析命令行参数
    args = parser.parse_args(('step1_splitter', '--input', ','.join(pic_list)))  # 测试
    # args = parser.parse_args()
    arg, _ = args.action

    # 运行第一步
    if arg == 'step1_splitter':
        for pic in args.input.split(','):
            pic = pic_list[4]
            image = MyImage(str(pic))
            image.get_con()
            image.split_by_con()
            # 保存中间结果
            image.draw_rect()  # 拆分的边框图
            image.draw_con()  # 识别到的pieces

            del image
            break

    # 运行第二步
    if arg == 'step2_concat':
        pass

    # 运行第三步
    # piece_list = [
    #     r'D:/workspace/MyUtils/fun/20211028-lcy-20X-1_Wholeslide_Default_Extended_image_results/pieces_results',
    #     r'D:/workspace/MyUtils/fun/20211028-lcy-20X-2_Wholeslide_Default_Extended_image_results/pieces_results',
    #     r'D:/workspace/MyUtils/fun/20211028-lcy-20X-3_Wholeslide_Default_Extended_image_results/pieces_results',
    #     r'D:/workspace/MyUtils/fun/20211028-lcy-20X-4_Wholeslide_Default_Extended_image_results/pieces_results',
    #     r'D:/workspace/MyUtils/fun/20211028-lcy-20X-5_Wholeslide_Default_Extended_image_results/pieces_results']
    #
    # MyImageTools.adjust_piece(*piece_list)
