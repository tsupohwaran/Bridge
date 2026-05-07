clear all
global proj_path = "/Users/pohwaran/Doctorate/Paper/Bridge"
do "$proj_path/code/00_setup/data_paths.do"

* 这段代码的目标是构造企业层面的 exposure 指标。
* 直观上，exposure 衡量的是：某个企业在“通桥”情形下，
* 与黄岛区/即墨区 2013 年新增企业活动的空间接触强度有多大。
*
* 具体做法分成两步：
* 1. 先在镇级层面计算 2013 年相对于 2012 年的新增企业数量 new_entry_firm_num。
* 2. 再把企业到各镇的通勤时间 dzj_prime 与这个新增企业数量合并，
*    用 new_entry_firm_num 作为权重，对同一企业在同一区内所有镇的 dzj_prime 做加权平均。
*    得到的县区层面加权平均通勤时间，就是 exposure。
*
* 因而，exposure 的计算公式可以写成：
* exposure_{i,c} = sum_j(dzj_prime_{i,j} * new_entry_firm_num_j) / sum_j(new_entry_firm_num_j)
* 其中 i 表示企业，c 表示县区（黄岛区或即墨区），j 表示县区内的镇。
*
* 第一步：计算黄岛区和即墨区各镇在 2013 年的新增企业数。
use "$busin_data_path/BusinRegis_QingDao_13.dta", clear
append using "$busin_data_path/BusinRegis_QingDao_12.dta"

* 统一经纬度变量名，便于后续空间匹配。
ren (经度 纬度) (longitude latitude)

* 将企业注册点按经纬度落到镇级多边形中，识别每家企业所在的镇。
geoinpoly latitude longitude using "$geo_raw_path/shapefiles/town_coord.dta"
ren _ID ID
merge m:1 ID using "$geo_raw_path/shapefiles/town_db.dta", nogen
ren (乡 县) (town county)
drop ID 省 市 treat geom

* 只保留研究关注的两个区：黄岛区和即墨区。
keep if county == "黄岛区" | county == "即墨区"

* 在“镇-年份”层面统计企业数量。
bys town year: gen firm_num = _N
duplicates drop town year, force
drop if missing(town)
keep town county firm_num year

* 宽表化后，让每个镇同时拥有 2012 和 2013 两年的企业数。
reshape wide firm_num, i(town) j(year)

* 新增企业数 = 2013 年企业数 - 2012 年企业数。
* 这里把它理解为该镇新增经济活动强度，后面会作为 exposure 的权重。
gen new_entry_firm_num = firm_num2013 - firm_num2012

* 去掉无法计算新增量或新增量为负的镇。
drop if missing(new_entry_firm_num) | new_entry_firm_num < 0
save "$regression_temp_path/new_entry_firm_2013.dta", replace

* 第二步：把企业到镇的通勤时间，与镇级新增企业数合并起来，计算 exposure。
*
* commute_time_qingdao_07_20_python.csv 中，每一行表示：
* 某家企业 id 到某个镇 town 的通勤时间 dzj_prime（单位：分钟）。
* 这里使用的是通桥情形下的 travel_time_min。
import delimited "$geo_processed_path/commute_time_qingdao_07_20_python.csv", clear
keep firm_id town_id travel_time_min
ren (firm_id town_id travel_time_min) (id town dzj_prime)

* 将每个镇对应的新增企业数和所属县区 merge 进来。
* merge 之后，一行数据可以理解为：
* “企业 i 到镇 j 的通勤时间” + “镇 j 属于哪个区” + “镇 j 的新增企业数量是多少”。
merge m:1 town using "$regression_temp_path/new_entry_firm_2013.dta", keep(3) nogen

sort id town

* exposure 的核心计算在这里：
* 对每个企业 id、每个县区 county，取该企业到该县区所有镇的 dzj_prime，
* 并用各镇的新增企业数 new_entry_firm_num 作为权重做加权平均。
*
* 这样做的含义是：
* 1. 如果某个镇新增企业更多，这个镇在 exposure 中的权重就更高；
* 2. 因而 exposure 更强调企业对“新增进入活动更集中”的镇的可达性；
* 3. 最终 exposure 越小，表示企业到这些新增企业更活跃区域的通勤时间越短。
bys id county: egen exposure = wtmean(dzj_prime), weight(new_entry_firm_num)

* 现在同一个企业-县区组内的每一行 exposure 都相同，
* 所以去重后只保留企业在该县区的一条 exposure 记录。
duplicates drop id county, force
keep id county exposure

* 把长表转成宽表，使每家企业同时拥有两个 exposure 指标：
* 一个对应黄岛区，一个对应即墨区。
reshape wide exposure, i(id) j(county, string)
ren (exposure黄岛区 exposure即墨区) (exposure_huangdao exposure_jimo)

* 最终输出是企业层面的 exposure 数据：
* exposure_huangdao = 企业对黄岛区新增企业活动的暴露度
* exposure_jimo     = 企业对即墨区新增企业活动的暴露度
replace exposure_huangdao = 100 / exposure_huangdao
replace exposure_jimo = 100 / exposure_jimo
save "$regression_temp_path/firm_exposure_qingdao.dta", replace
