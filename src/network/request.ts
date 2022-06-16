import axios from "axios";
export function ins() {
  // 进行基础配置
  const instance = axios.create({
    // 在进行完跨域配置后，访问DRF提供的后端接口
    // 使用axios.create进行基本配置
    baseURL: "http://127.0.0.1:8000/",
    timeout: 10000,
  });
  return instance;
}
